package nl.systemsgenetics.downstreamer.runners;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import htsjdk.samtools.util.IntervalTreeMap;
import me.tongfei.progressbar.ProgressBar;
import nl.systemsgenetics.downstreamer.DownstreamerOptions;
import nl.systemsgenetics.downstreamer.DownstreamerStep2Results;
import nl.systemsgenetics.downstreamer.DownstreamerStep3Results;
import nl.systemsgenetics.downstreamer.gene.Gene;
import nl.systemsgenetics.downstreamer.io.ExcelWriter;
import nl.systemsgenetics.downstreamer.io.IoUtils;
import nl.systemsgenetics.downstreamer.pathway.PathwayDatabase;
import nl.systemsgenetics.downstreamer.pathway.PathwayEnrichments;
import nl.systemsgenetics.downstreamer.summarystatistic.LdScore;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.descriptive.moment.Skewness;
import org.apache.commons.math3.stat.inference.TTest;
import org.apache.log4j.Logger;
import org.molgenis.genotype.util.Ld;
import umcg.genetica.graphics.panels.HistogramPanel;
import umcg.genetica.math.PcaColt;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.PearsonRToZscoreBinned;
import umcg.genetica.math.stats.ZScores;

import java.io.*;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.IntStream;
import java.util.zip.GZIPInputStream;

import nl.systemsgenetics.downstreamer.containers.LeadVariant;
import umcg.genetica.collections.ChrPosTreeMap;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetFastSubsetLoader;

/**
 * Collection of runners that handle pre-processing steps for Depict 2 analysis.
 */
public class DownstreamerUtilities {

    private static final Logger LOGGER = Logger.getLogger(DownstreamerUtilities.class);

    private static HashMap<String, HashMap<String, NearestVariant>> traitGeneDist = null;

    /**
     * Create a gene gene correlation matrix based on a (eigenvector) matrix.
     *
     * @param options
     * @throws Exception
     */
    public static void correlateGenes(DownstreamerOptions options) throws FileNotFoundException, Exception {

        DoubleMatrixDataset<String, String> expressionMatrix;

        if (options.getGeneInfoFile() != null) {
            final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
            final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(options.getGeneInfoFile()))).withCSVParser(parser).withSkipLines(0).build();

            final HashSet<String> genes = new HashSet<>();
            String[] nextLine;
            while ((nextLine = reader.readNext()) != null) {
                genes.add(nextLine[0]);
            }
            LOGGER.info("Read " + genes.size() + " genes to load");
            if (options.getGwasZscoreMatrixPath().endsWith(".txt") || options.getGwasZscoreMatrixPath().endsWith("txt.gz")) {
                expressionMatrix = DoubleMatrixDataset.loadSubsetOfTextDoubleData(options.getGwasZscoreMatrixPath(), '\t', genes, null);
            } else {

                DoubleMatrixDatasetFastSubsetLoader loader = new DoubleMatrixDatasetFastSubsetLoader(options.getGwasZscoreMatrixPath());
                Map<String, Integer> rows = loader.getOriginalRowMap();

                genes.retainAll(rows.keySet());

                expressionMatrix = loader.loadSubsetOfRowsBinaryDoubleData(genes);
            }

        } else {
            if (options.getGwasZscoreMatrixPath().endsWith(".txt") || options.getGwasZscoreMatrixPath().endsWith("txt.gz")) {
                expressionMatrix = DoubleMatrixDataset.loadDoubleTextData(options.getGwasZscoreMatrixPath(), '\t');
            } else {
                expressionMatrix = DoubleMatrixDataset.loadDoubleBinaryData(options.getGwasZscoreMatrixPath());
            }
        }

        // Normalize the input data
        if (options.isNormalizeEigenvectors()) {
            expressionMatrix.normalizeRows();
            expressionMatrix.normalizeColumns();
            LOGGER.info("Data row normalized and then column normalized");
        }

        // Optionally select a subset of columns to use
        String[] cols = options.getColumnsToExtract();
        if (cols != null) {
            Set<String> columnsToExtract = new HashSet<>();
            for (String colname : cols) {
                if (expressionMatrix.getColObjects().contains(colname)) {
                    columnsToExtract.add(colname);
                } else {
                    LOGGER.warn(colname + " is missing in input matrix, ommiting col in output");
                }
            }

            expressionMatrix = expressionMatrix.viewColSelection(columnsToExtract);
        }

        // Calculate the correlation matrix
        LOGGER.info("Loaded expression matrix with " + expressionMatrix.rows() + " genes and " + expressionMatrix.columns() + " observations");
        DoubleMatrixDataset<String, String> corMatrix = expressionMatrix.viewDice().calculateCorrelationMatrix();
        LOGGER.info("Done calculating correlations");

        // Convert Pearson R to Z-scores
        if (options.isCorMatrixZscores()) {
            PearsonRToZscoreBinned r2zScore = new PearsonRToZscoreBinned(10000000, expressionMatrix.columns());
            r2zScore.inplaceRToZ(corMatrix);
            LOGGER.info("Converted correlations to Z-scores");
        }

        // Set diagonal of matrix to zero
        for (int i = 0; i < corMatrix.columns(); ++i) {
            corMatrix.setElementQuick(i, i, 0);
        }
        LOGGER.info("Diagnonal set to zero as this might inflate coregulation towards genes in GWAS loci");

        // Save
        corMatrix.saveBinary(options.getOutputBasePath());
        LOGGER.info("Correlation matrix saved to: " + options.getOutputBasePath() + ".dat");

        // Calculate per gene distribution metrics
        LOGGER.info("Calculating per gene distribution metrics");
        DoubleMatrixDataset<String, String> perGeneDistMetrics = DownstreamerUtilities.calculateDistributionMetricsPerRow(corMatrix);
        perGeneDistMetrics.save(options.getOutputBasePath() + ".coregulation.dist.metrics.txt");
        LOGGER.info("Done");

    }

    /**
     * Converts a matrix of Pearson R values to z-scores. Diagonal of this
     * matrix is set to zero.
     *
     * @param options
     * @throws FileNotFoundException
     * @throws Exception
     */
    public static void convertRtoZscore(DownstreamerOptions options) throws FileNotFoundException, Exception {
        DoubleMatrixDataset<String, String> corMatrix;

        if (options.getGwasZscoreMatrixPath().endsWith(".txt") || options.getGwasZscoreMatrixPath().endsWith("txt.gz")) {
            corMatrix = DoubleMatrixDataset.loadDoubleTextData(options.getGwasZscoreMatrixPath(), '\t');
        } else {
            corMatrix = DoubleMatrixDataset.loadDoubleBinaryData(options.getGwasZscoreMatrixPath());
        }

        PearsonRToZscoreBinned r2zScore = new PearsonRToZscoreBinned(10000000, options.getNumberSamplesUsedForCor());
        r2zScore.inplaceRToZ(corMatrix);
        LOGGER.info("Converted correlations to Z-scores");

        for (int i = 0; i < corMatrix.columns(); ++i) {
            corMatrix.setElementQuick(i, i, 0);
        }
        LOGGER.info("Diagnonal set to zero as this might inflate coregulation towards genes in GWAS loci");

        corMatrix.saveBinary(options.getOutputBasePath());

        LOGGER.info("Correlation matrix saved to: " + options.getOutputBasePath() + ".dat");

    }

    /**
     * Run PCA analysis on a binary matrix using PcaColt.
     *
     * @param options
     * @throws IOException
     */
    public static void doPcaOnBinMatrix(DownstreamerOptions options) throws IOException {

        final DoubleMatrixDataset<String, String> dataset = DoubleMatrixDataset.loadDoubleBinaryData(options.getGwasZscoreMatrixPath());

        //if debug is enabled keep cov matrix in memory
        PcaColt pcaRes = new PcaColt(dataset, true, true, LOGGER.isDebugEnabled());

        pcaRes.getEigenvectors().save(options.getOutputBasePath() + "_eigenVectors.txt");
        pcaRes.getEigenValues().save(options.getOutputBasePath() + "_eigenValues.txt");
        pcaRes.getPcs().save(options.getOutputBasePath() + "_pcs.txt");
        if (LOGGER.isDebugEnabled()) {
            pcaRes.getCovMatrix().save(options.getOutputBasePath() + "_correlationMatrix.txt");
        }

    }

    /**
     * Utility to get the force normalized gene pvalues with tie resolving.
     * Shares some duplicate code with Depict2MainAnalysis.run2(). This can be
     * cleaned up in future
     *
     * @param options
     * @throws Exception
     */
    public static void getNormalizedGwasGenePvalues(DownstreamerOptions options) throws Exception {
        getNormalizedGwasGenePvaluesReturn(options);
    }

    /**
     * Utility to get the force normalized gene pvalues with tie resolving.
     * Shares some duplicate code with Depict2MainAnalysis.run2(). This can be
     * cleaned up in future
     *
     * @param options
     * @throws Exception
     */
    public static DoubleMatrixDataset<String, String> getNormalizedGwasGenePvaluesReturn(DownstreamerOptions options) throws Exception {

        DoubleMatrixDataset<String, String> genePvalues;
        List<Gene> genes;
        DoubleMatrixDataset<String, String> geneVariantCount;
        DoubleMatrixDataset<String, String> geneMaxSnpZscore;

        LOGGER.info("Continuing previous analysis by loading gene p-values");
        if (new File(options.getRun1BasePath() + "_genePvalues.dat").exists()) {
            genePvalues = DoubleMatrixDataset.loadDoubleBinaryData(options.getRun1BasePath() + "_genePvalues");
            // Always load to avoid nullpointers
            geneMaxSnpZscore = DoubleMatrixDataset.loadDoubleBinaryData(options.getRun1BasePath() + "_geneMaxSnpScores");

        } else {
            LOGGER.fatal("Could not find gene pvalues at: " + options.getRun1BasePath() + "_genePvalues.dat");
            LOGGER.fatal("First use --mode RUN to calculate gene p-values");
            return null;
        }

        geneVariantCount = DoubleMatrixDataset.loadDoubleTextData(options.getRun1BasePath() + "_geneVariantCount.txt", '\t');
        LOGGER.info("Gene p-values loaded");
        genes = IoUtils.readGenes(options.getGeneInfoFile());
        LOGGER.info("Loaded " + genes.size() + " genes");

        // Identify genes with at least one variant in window
        final HashSet<String> selectedGenes = new HashSet<>();
        final ArrayList<String> allGenes = geneVariantCount.getRowObjects();
        final int totalGeneCount = allGenes.size();
        for (int g = 0; g < totalGeneCount; ++g) {
            if (geneVariantCount.getElementQuick(g, 0) > 0) {
                selectedGenes.add(allGenes.get(g));
            }
        }

        // Select genes that have a gene pvalue, to avoid issues with the normalization, and to keep consistency
        // with the PathwayEnrichments.
        genePvalues = genePvalues.viewRowSelection(selectedGenes);

        LOGGER.info(genePvalues.rows() + " have a gene pvalue");
        final DoubleMatrix2D matrix = genePvalues.getMatrix();

        // Inplace convert gene p-values to z-scores
        IntStream.range(0, matrix.rows()).parallel().forEach(r -> {
            for (int c = 0; c < matrix.columns(); ++c) {
                matrix.setQuick(r, c, -ZScores.pToZTwoTailed(matrix.getQuick(r, c)));
            }
        });

        LOGGER.info("Force normalizing gene p-values / z-scores");
        DoubleMatrixDataset<String, String> normalizedGwasGeneScores;
        normalizedGwasGeneScores = PathwayEnrichments.createColumnForceNormalDuplicate(genePvalues, geneMaxSnpZscore);
        normalizedGwasGeneScores.save(options.getOutputBasePath() + "_normalizedGenePvalues.txt");

        return (normalizedGwasGeneScores);
    }

    /**
     * Re-generate the excel file from existing files in the intermediate
     * folder.
     *
     * @param options
     * @throws Exception
     */
    public static void generateExcelFromIntermediates(DownstreamerOptions options) throws Exception {

        DownstreamerStep2Results step2 = loadExistingStep2Results(options);
        ExcelWriter writer = new ExcelWriter(step2.getGenePvalues().getColObjects(), options);

        writer.saveStep2Excel(step2);
        //writer.saveGenePvalueExcel(step2.getGenePvalues());

        if (options.getPathwayDatabasesToAnnotateWithGwas().size() >= 1) {
            DownstreamerStep3Results step3 = DownstreamerMainAnalysis.step3(options);
            writer.saveStep3Excel(step2, step3);
        }

    }

    /**
     * Generate an excel file with the z-scores of the pathways for all bonf.
     * sig. genes and pathways.
     *
     * @param options
     * @throws Exception
     */
    public static void generatePathwayLoadingExcel(DownstreamerOptions options) throws Exception {
        DownstreamerStep2Results step2 = loadExistingStep2Results(options);

        ExcelWriter writer = new ExcelWriter(step2.getGenePvalues().getColObjects(), options);
        writer.savePathwayLoadings(step2);
    }

    /**
     * Load existing results from step 2 from storage
     *
     * @param options
     * @return
     * @throws Exception
     */
    public static DownstreamerStep2Results loadExistingStep2Results(DownstreamerOptions options) throws Exception {
        return loadExistingStep2Results(options, false);
    }

    /**
     * Load existing results from step 2 from storage. If matchToNormalizedPvalues = true, the genePvalues and
     * normalizedGenePvalues are matched.
     *
     * @param options
     * @return
     * @throws Exception
     */
    public static DownstreamerStep2Results loadExistingStep2Results(DownstreamerOptions options, boolean matchToNormalizedPvalues) throws Exception {

        DoubleMatrixDataset<String, String> genePvalues = DoubleMatrixDataset.loadDoubleBinaryData(options.getRun1BasePath() + "_genePvalues");
        DoubleMatrixDataset<String, String> normalizedGenePvalues;
        if (options.isForceNormalGenePvalues()) {
            normalizedGenePvalues = getNormalizedGwasGenePvaluesReturn(options);
        } else {
            normalizedGenePvalues = null;
        }

        if (matchToNormalizedPvalues) {
            genePvalues = genePvalues.viewRowSelection(normalizedGenePvalues.getRowObjects());
        }

        final List<PathwayDatabase> pathwayDatabases = options.getPathwayDatabases();
        ArrayList<PathwayEnrichments> pathwayEnrichments = new ArrayList<>(pathwayDatabases.size());
        for (PathwayDatabase pathwayDatabase : pathwayDatabases) {
            pathwayEnrichments.add(new PathwayEnrichments(pathwayDatabase, options.getIntermediateFolder(), options.isExcludeHla()));
        }
        return new DownstreamerStep2Results(pathwayEnrichments, genePvalues, normalizedGenePvalues);

    }

    /**
     * Load a co-regulation matrix and set any gene-gene correlation between genes closer than 250kb to zero.
     *
     * @param options
     * @throws IOException
     * @throws Exception
     */
    public static void removeLocalGeneCorrelations(DownstreamerOptions options) throws IOException, Exception {

        LinkedHashMap<String, Gene> genes = IoUtils.readGenesMap(options.getGeneInfoFile());
        LOGGER.info("Loaded " + genes.size() + " genes");

        DoubleMatrixDataset<String, String> corMatrix = DoubleMatrixDataset.loadDoubleBinaryData(options.getGwasZscoreMatrixPath());

        List<String> genesToKeep = corMatrix.getRowObjects();
        LOGGER.info("Read " + genesToKeep.size() + " genes in correlation matrix");
        genesToKeep.retainAll(genes.keySet());
        LOGGER.info("Retained " + genesToKeep.size() + " genes that overlap with --genes file");
        corMatrix = corMatrix.viewSelection(genesToKeep, genesToKeep);


        if (!corMatrix.getHashRows().keySet().containsAll(corMatrix.getHashCols().keySet())) {
            throw new Exception("Co-expression matrix is not squared with same row and col names");
        }

        if (!genes.keySet().containsAll(corMatrix.getHashRows().keySet())) {
            throw new Exception("Not all genes Co-expression matrix are found in gene mapping file");
        }

        final int genesInMatrix = corMatrix.rows();
        final ArrayList<String> geneOrder = corMatrix.getRowObjects();

        int overlappingGenePairs = 0;

        for (int i = 0; i < genesInMatrix; ++i) {

            //diagnoal always 0
            corMatrix.setElementQuick(i, i, 0);

            Gene geneI = genes.get(geneOrder.get(i));

            for (int j = i + 1; j < genesInMatrix; ++j) {

                Gene geneJ = genes.get(geneOrder.get(j));

                if (geneI.withinDistanceOf(geneJ, options.getCisWindowExtend())) {
                    corMatrix.setElementQuick(i, j, 0);
                    corMatrix.setElementQuick(j, i, 0);
                    ++overlappingGenePairs;
                }

            }

        }

        LOGGER.info("Identified " + overlappingGenePairs + " overlapping gene-gene pairs within " + options.getCisWindowExtend() + "b.");

        corMatrix.saveBinary(options.getOutputBasePath());

    }


    /**
     * Calculate skewness, kurtosis, mean and variance of the null distribution for a pathway database
     *
     * @throws Exception
     */
    public static DoubleMatrixDataset<String, String> calculateDistributionMetricsPerRow(DoubleMatrixDataset<String, String> data) {

        List<String> colnames = new ArrayList<>();
        colnames.add("mean");
        colnames.add("sd");
        colnames.add("skewness");
        colnames.add("kurtosis");
        colnames.add("max");
        colnames.add("min");

        DoubleMatrixDataset<String, String> outputStats = new DoubleMatrixDataset<>(data.getRowObjects(), colnames);

        for (int r = 0; r < data.rows(); r++) {
            DescriptiveStatistics curDesc = new DescriptiveStatistics(data.getRow(r).toArray());

            int c = 0;
            outputStats.setElementQuick(r, c++, curDesc.getMean());
            outputStats.setElementQuick(r, c++, curDesc.getStandardDeviation());
            outputStats.setElementQuick(r, c++, curDesc.getSkewness());
            outputStats.setElementQuick(r, c++, curDesc.getKurtosis());
            outputStats.setElementQuick(r, c++, curDesc.getMax());
            outputStats.setElementQuick(r, c++, curDesc.getMin());
        }

        return outputStats;
    }


    public static void calculateMeanLdScorePerGene(DownstreamerOptions options) throws IOException {

        List<Gene> genes = IoUtils.readGenes(options.getGeneInfoFile());
        int window = options.getWindowExtend();

        IntervalTreeMap<LdScore> ldScores = readLdScores(options);


        for (Gene curGene : genes) {


        }


    }

    /**
     * Input a normalized expression matrix, and perform a t-test between all the samples indicated in the grouping file and the rest.
     * Is parallelizable.
     *
     * @param options
     */
    public static void generateMarkerGenes(DownstreamerOptions options) throws Exception {

        Set<String> allListedSamples = new HashSet<>();
        Map<String, Set<String>> sampleGroups = new HashMap<>();
        BufferedReader reader = new BufferedReader(new FileReader(options.getX()));
/*		while (reader.ready()) {
			String[] line = reader.readLine().split("\t");

			Set<String> curSamples = new HashSet<>();
			Collections.addAll(curSamples, Arrays.copyOfRange(line, 1, line.length));

			sampleGroups.put(line[0], curSamples);
			allListedSamples.addAll(curSamples);
		}*/

        Set<String> availIds = new HashSet<>(DoubleMatrixDataset.readDoubleTextDataRowNames(options.getGwasZscoreMatrixPath() + ".cols.txt", '\t'));
        LOGGER.info("Read " + availIds.size() + " available samples");

        while (reader.ready()) {
            String[] line = reader.readLine().split("\t");
            if (line.length != 2) {
                throw new IllegalArgumentException("Line in grouping file contains != 2 columns");
            }

            if (availIds.contains(line[1])) {
                if (sampleGroups.containsKey(line[0])) {
                    sampleGroups.get(line[0]).add(line[1]);
                } else {
                    Set<String> curSet = new HashSet<>();
                    curSet.add(line[1]);
                    sampleGroups.put(line[0], curSet);
                }
                allListedSamples.add(line[1]);
            }

        }


        LOGGER.info("Read " + sampleGroups.size() + " sample groups over " + allListedSamples.size() + " samples");

        List<Gene> genes = IoUtils.readGenes(options.getGeneInfoFile());
        Set<String> geneIds = new HashSet<>();
        for (Gene gene : genes) {
            geneIds.add(gene.getGene());
        }

        LOGGER.info("Read " + genes.size() + " genes groups to test");

        // Bulk load expression matrix into memory, could convert to streaming algorithm for more mem efficiency
        LOGGER.info("Loading expression data");

        geneIds.retainAll(DoubleMatrixDataset.readDoubleTextDataRowNames(options.getGwasZscoreMatrixPath() + ".rows.txt", '\t'));

        final DoubleMatrixDataset<String, String> expression = DoubleMatrixDataset.loadSubsetOfRowsBinaryDoubleData(options.getGwasZscoreMatrixPath(), geneIds).viewColSelection(allListedSamples);
        LOGGER.info("Done Loading expression data");

        // Ttest object, and storage for output
        final TTest tTest = new TTest();
        final DoubleMatrixDataset<String, String> pvalues = new DoubleMatrixDataset<>(geneIds, sampleGroups.keySet());
        final DoubleMatrixDataset<String, String> tstats = new DoubleMatrixDataset<>(geneIds, sampleGroups.keySet());


        LOGGER.info("Running T-tests");
        final ProgressBar bp = new ProgressBar("genes", geneIds.size());
        bp.step();
        // Loop over genes, to make more memory efficient, could load one gene at the time
        // Currently with expression matrix of 19x17k its quite doable
        //for (String gene : geneIds) {

        ForkJoinPool forkJoinPool = null;
        try {
            forkJoinPool = new ForkJoinPool(options.getNumberOfThreadsToUse());
            forkJoinPool.submit(() -> geneIds.parallelStream().forEach(gene -> {
                Set<String> curGene = new HashSet<>();
                curGene.add(gene);

                // Loop over the differnet tissue groups
                for (String tissue : sampleGroups.keySet()) {

                    // Determine the groups to test
                    Set<String> groupA = sampleGroups.get(tissue);
                    Set<String> groupB = new HashSet<>();

                    for (String tissue2 : sampleGroups.keySet()) {
                        if (!tissue2.equals(tissue)) {
                            groupB.addAll(sampleGroups.get(tissue2));
                        }
                    }

                    // Subset the 2 sets
                    DoubleMatrix1D groupAValues = expression.viewSelection(curGene, groupA).getRow(gene);
                    DoubleMatrix1D groupBValues = expression.viewSelection(curGene, groupB).getRow(gene);

                    // Run T-test
                    double tstat = tTest.t(groupAValues.toArray(), groupBValues.toArray());
                    double pval = tTest.tTest(groupAValues.toArray(), groupBValues.toArray());

                    // Save to output matrix
                    pvalues.setElement(gene, tissue, pval);
                    tstats.setElement(gene, tissue, tstat);
                }

                bp.step();
            })).get();

        } catch (Exception e) {
            throw new RuntimeException(e);
        } finally {
            if (forkJoinPool != null) {
                forkJoinPool.shutdown();
                bp.close();
                LOGGER.info("Done calculating T-tests");
            }
        }

        // Save output
        pvalues.save(options.getOutputBasePath() + ".pvalues.txt");
        tstats.save(options.getOutputBasePath() + ".tstats.txt");
        LOGGER.info("Done");

    }

    /**
     * HashMap<Trait, Map<Gene, distance>
     * <0 for other chr or outside cis window
     *
     * @param options
     * @return
     * @throws IOException
     */
    public static HashMap<String, HashMap<String, NearestVariant>> getDistanceGeneToTopCisSnpPerTrait(
            final DownstreamerOptions options) throws IOException {

        if (traitGeneDist != null) {
            return traitGeneDist;
        }


        Map<String, ChrPosTreeMap<LeadVariant>> indepVariantsAsSummaryStatisticsRecord = IoUtils.loadLeadVariantsPerTrait(options);
        LinkedHashMap<String, Gene> genes = IoUtils.readGenesMap(options.getGeneInfoFile());
        final int cisExtent = options.getCisWindowExtend();

        HashMap<String, HashMap<String, NearestVariant>> traitGeneDist2 = new HashMap<>(indepVariantsAsSummaryStatisticsRecord.size());

        for (Map.Entry<String, ChrPosTreeMap<LeadVariant>> traitEntry : indepVariantsAsSummaryStatisticsRecord.entrySet()) {

            String trait = traitEntry.getKey();
            ChrPosTreeMap<LeadVariant> topHits = traitEntry.getValue();

            HashMap<String, NearestVariant> geneDist = new HashMap<>(genes.size());
            traitGeneDist2.put(trait, geneDist);

            for (Gene gene : genes.values()) {

                int geneStart = Math.min(gene.getStart(), gene.getEnd());
                int geneEnd = Math.max(gene.getStart(), gene.getEnd());

                int minDist = Integer.MAX_VALUE;
                LeadVariant nearestVariant = null;

                for (LeadVariant cisVariant : topHits.getChrRange(gene.getContig(), geneStart - cisExtent, true, geneEnd + cisExtent, true).values()) {
                    if (cisVariant.getPos() >= geneStart && cisVariant.getPos() <= geneEnd) {
                        minDist = 0;
                        nearestVariant = cisVariant;
                        continue;
                    }

                    int dist;
                    if (cisVariant.getPos() < geneStart) {
                        dist = geneStart - cisVariant.getPos();
                    } else {
                        dist = cisVariant.getPos() - geneEnd;
                    }
                    if (dist < minDist) {
                        minDist = dist;
                        nearestVariant = cisVariant;
                    }

                }

                if (minDist == Integer.MAX_VALUE) {
                    //no variant found in window
                    minDist = -9;
                }

                geneDist.put(gene.getGene(), new NearestVariant(nearestVariant, minDist));

            }


        }

        traitGeneDist = traitGeneDist2;

        return traitGeneDist2;

    }

    public static class NearestVariant {

        private final LeadVariant nearestVariant;
        private final int distance;

        public NearestVariant(LeadVariant nearestVariant, int distance) {
            this.nearestVariant = nearestVariant;
            this.distance = distance;
        }

        public LeadVariant getNearestVariant() {
            return nearestVariant;
        }

        public int getDistance() {
            return distance;
        }

    }

    /**
     * Read LD score files into an IntervalTreeMap in the format provided by https://github.com/bulik/ldsc
     * @param options
     * @return
     * @throws IOException
     */
    public static IntervalTreeMap<LdScore> readLdScores(DownstreamerOptions options) throws IOException {

        IntervalTreeMap<LdScore> output = new IntervalTreeMap<>();
        for (int i=1; i < 23; i++) {
            BufferedReader reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(options.getGwasZscoreMatrixPath() + "/" + i + ".l2.ldscore.gz"))));
            String[] line = reader.readLine().split("\t");

            if (line[0].equals("CHR")) {
                continue;
            }

            LdScore curLdscore =  new LdScore(line[0],
                    Integer.parseInt(line[2]),
                    line[1],
                    Double.parseDouble(line[5]));

            output.put(curLdscore, curLdscore);
        }

        return output;
    }

    /**
     * Calculate the benjamini hochberg adjusted p-values form a doublematrixdataset. Preserves the order of the orginal
     * input.
     * Adapted from: https://github.com/cBioPortal/cbioportal/blob/master/core/src/main/java/org/mskcc/cbio/portal/stats/BenjaminiHochbergFDR.java
     * and
     * https://stats.stackexchange.com/questions/238458/whats-the-formula-for-the-benjamini-hochberg-adjusted-p-value
     *
     * @param pvalues
     * @return
     */
    public static DoubleMatrixDataset<String, String> adjustPvaluesBenjaminiHochberg(DoubleMatrixDataset<String, String> rawPvalues) {
        DoubleMatrixDataset<String, String> output = new DoubleMatrixDataset<>(rawPvalues.getRowObjects(), rawPvalues.getColObjects());
        int m = rawPvalues.rows();

        for (int c = 0; c < rawPvalues.columns(); c++) {

            String colname = output.getColObjects().get(c);

            // Sort the p-values preserving the ids
            List<DoubleElement> sortedPvalues = new ArrayList<>(rawPvalues.rows());
            for (int r = 0; r < rawPvalues.rows(); r++) {
                sortedPvalues.add(new DoubleElement(rawPvalues.getElementQuick(r, c), rawPvalues.getRowObjects().get(r)));
            }
            sortedPvalues.sort(Comparator.comparing(DoubleElement::getValue));

            List<DoubleElement> adjustedPvalues = new ArrayList<>(sortedPvalues);

            // iterate through all p-values:  largest to smallest
            for (int i = m - 1; i >= 0; i--) {
                if (i == m - 1) {
                    adjustedPvalues.set(i, sortedPvalues.get(i));
                } else {
                    double unadjustedPvalue = sortedPvalues.get(i).value;
                    int divideByM = i + 1;
                    double left = adjustedPvalues.get(i + 1).value;
                    double right = (m / (double) divideByM) * unadjustedPvalue;
                    adjustedPvalues.set(i, new DoubleElement(Math.min(left, right), sortedPvalues.get(i).id));
                }
            }

            for (DoubleElement curElement : adjustedPvalues) {
                output.setElement(curElement.id, colname, curElement.value);
            }

        }

        return output;
    }


    /**
     * Links a double with an ID. Used in adjustPvaluesBenjaminiHochberg.
     */
    private static class DoubleElement {
        protected double value;
        protected String id;

        public DoubleElement(double value, String id) {
            this.value = value;
            this.id = id;
        }

        public double getValue() {
            return value;
        }

        public String getId() {
            return id;
        }
    }

}
