package nl.systemsgenetics.downstreamer.runners;

import nl.systemsgenetics.downstreamer.DownstreamerStep3Results;
import nl.systemsgenetics.downstreamer.DownstreamerStep2Results;
import nl.systemsgenetics.downstreamer.DownstreamerMode;
import nl.systemsgenetics.downstreamer.DownstreamerStep1Results;
import nl.systemsgenetics.downstreamer.Downstreamer;
import nl.systemsgenetics.downstreamer.DownstreamerOptions;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import nl.systemsgenetics.downstreamer.gene.Gene;
import nl.systemsgenetics.downstreamer.gene.GenePvalueCalculator;
import nl.systemsgenetics.downstreamer.io.IoUtils;
import nl.systemsgenetics.downstreamer.pathway.PathwayDatabase;
import nl.systemsgenetics.downstreamer.pathway.PathwayEnrichments;
import nl.systemsgenetics.downstreamer.summarystatistic.Locus;
import nl.systemsgenetics.downstreamer.summarystatistic.LocusUtils;
import nl.systemsgenetics.downstreamer.summarystatistic.SummaryStatisticRecord;
import org.apache.log4j.Logger;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.ZScores;

import java.io.*;
import java.util.*;
import java.util.concurrent.ThreadLocalRandom;
import java.util.stream.IntStream;

/**
 * Runners for the main Depict 2 analysis.
 */
public class DownstreamerMainAnalysis {

    private static final Logger LOGGER = Logger.getLogger(DownstreamerMainAnalysis.class);

    /**
     * Calculate the real gene p-values from GWAS z-scores as well as those for x ammount of random GWAS phenotypes.
     *
     * @param options
     * @return DownstreamerStep1Results
     * @throws IOException
     * @throws Exception
     */
    public static DownstreamerStep1Results step1(DownstreamerOptions options) throws IOException, Exception {

        //Test here to prevent chrash after running for a while
        if (!new File(options.getGwasZscoreMatrixPath() + ".dat").exists()) {
            throw new FileNotFoundException("GWAS matrix does not exist at: " + options.getGwasZscoreMatrixPath() + ".dat");
        }

        RandomAccessGenotypeData referenceGenotypeData = IoUtils.readReferenceGenotypeDataMatchingGwasSnps(options);
        LOGGER.info("Done loading genotype data");

        List<Gene> genes = IoUtils.readGenes(options.getGeneInfoFile());
        LOGGER.info("Loaded " + genes.size() + " genes");

        double[] randomChi2 = generateRandomChi2(options.getNumberOfPermutationsRescue(), 500);

        LOGGER.info("Prepared reference null distribution with " + Downstreamer.LARGE_INT_FORMAT.format(randomChi2.length) + " values");

        File usedVariantsPerGeneFile = options.isSaveUsedVariantsPerGene() ? new File(options.getOutputBasePath() + "_usedVariantsPerGene.txt") : null;

        GenePvalueCalculator gpc = new GenePvalueCalculator(options.getGwasZscoreMatrixPath(),
                referenceGenotypeData,
                genes,
                options.getWindowExtend(),
                options.getMaxRBetweenVariants(),
                options.getNumberOfPermutations(),
                options.getNumberOfPermutationsRescue(),
                options.getOutputBasePath(),
                randomChi2,
                options.correctForLambdaInflation(),
                options.getPermutationGeneCorrelations(),
                options.getPermutationPathwayEnrichment() + options.getPermutationFDR(),
                options.getDebugFolder(),
                options.getVariantGeneLinkingFile(),
                usedVariantsPerGeneFile);

        DoubleMatrixDataset<String, String> genePvalues = gpc.getGenePvalues();
        DoubleMatrixDataset<String, String> genePvaluesNullGwas = gpc.getGenePvaluesNullGwas();
        DoubleMatrixDataset<String, String> geneVariantCount = gpc.getGeneVariantCount();
        DoubleMatrixDataset<String, String> geneMaxSnpZscore = gpc.getGeneMaxSnpZscore();
        DoubleMatrixDataset<String, String> geneMaxSnpZscoreNullGwas = gpc.getGeneMaxSnpZscoreNullGwas();

        genePvalues.saveBinary(options.getOutputBasePath() + "_genePvalues");
        genePvaluesNullGwas.saveBinary(options.getOutputBasePath() + "_genePvaluesNullGwas");
        geneVariantCount.save(options.getOutputBasePath() + "_geneVariantCount.txt");
        geneMaxSnpZscore.saveBinary(options.getOutputBasePath() + "_geneMaxSnpScores");
        geneMaxSnpZscoreNullGwas.saveBinary(options.getOutputBasePath() + "_geneMaxSnpZscoresNullGwas");

        if (LOGGER.isDebugEnabled()) {
            gpc.getGeneMaxPermutationCount().save(options.getOutputBasePath() + "_geneMaxPermutationUsed.txt");
            gpc.getGeneRuntime().save(options.getOutputBasePath() + "_geneRuntime.txt");
        }
        LOGGER.info("Gene p-values saved. If needed the analysis can be resummed from this point using --mode STEP2 and exactly the same output path and genes file");

        PruneToIndependentTopHits.prune(options, referenceGenotypeData);

        return new DownstreamerStep1Results(genePvalues, genePvaluesNullGwas, geneVariantCount, geneMaxSnpZscore, geneMaxSnpZscoreNullGwas);
    }

    /**
     * Calculate the association between gene pvalues and pathway databases using a GLS model that accounts
     * for the correlation between gene-pvalues.
     *
     * @param options
     * @param step1Res
     * @return DownstreamerStep2Results
     * @throws IOException
     * @throws Exception
     */
    public static DownstreamerStep2Results step2(DownstreamerOptions options, DownstreamerStep1Results step1Res) throws IOException, Exception {

        options.getIntermediateFolder().mkdir();
		
		final List<PathwayDatabase> pathwayDatabases = options.getPathwayDatabases();
		
		for(PathwayDatabase pd : pathwayDatabases){
			if(!pd.exist()){
				throw new FileNotFoundException("Could not read: " + pd.getLocation() + ".dat");
			}
		}

        if (options.getMode() == DownstreamerMode.STEP2) {
            LOGGER.info("Continuing previous analysis by loading gene p-values");
            step1Res = DownstreamerStep1Results.loadFromDisk(options.getRun1BasePath());
            LOGGER.info("Gene p-values loaded");
        }

        DoubleMatrixDataset<String, String> genePvalues = step1Res.getGenePvalues();
        DoubleMatrixDataset<String, String> genePvaluesNullGwas = step1Res.getGenePvaluesNullGwas();
        DoubleMatrixDataset<String, String> geneVariantCount = step1Res.getGeneVariantCount();
        DoubleMatrixDataset<String, String> geneMaxSnpZscore = step1Res.getGeneMaxSnpZscore();
        DoubleMatrixDataset<String, String> geneMaxSnpZscoreNullGwas = step1Res.getGeneMaxSnpZscoreNullGwas();
        LinkedHashMap<String, Gene> genes = IoUtils.readGenesMap(options.getGeneInfoFile());

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

        final DoubleMatrix2D matrix = genePvalues.getMatrix();

        // Inplace convert gene p-values to z-scores
        IntStream.range(0, matrix.rows()).parallel().forEach(r -> {
            for (int c = 0; c < matrix.columns(); ++c) {
                matrix.setQuick(r, c, -ZScores.pToZTwoTailed(matrix.getQuick(r, c)));
            }
        });

        DoubleMatrix2D matrixNull = genePvaluesNullGwas.getMatrix();

        IntStream.range(0, matrixNull.rows()).parallel().forEach(r -> {
            for (int c = 0; c < matrixNull.columns(); ++c) {
                matrixNull.setQuick(r, c, -ZScores.pToZTwoTailed(matrixNull.getQuick(r, c)));
            }
        });

        LOGGER.info("Number of genes with atleast one variant in specified window: " + Downstreamer.LARGE_INT_FORMAT.format(selectedGenes.size()));
        final HashSet<String> hlaGenes;
        if (options.isExcludeHla()) {
            hlaGenes = new HashSet<>();
            for (Gene gene : genes.values()) {
                if (gene.getChr().equals("6") && ((gene.getStart() > 20000000 && gene.getStart() < 40000000) || (gene.getEnd() > 20000000 && gene.getEnd() < 40000000))) {
                    hlaGenes.add(gene.getGene());
                }
            }
            LOGGER.info("Excluding " + hlaGenes.size() + " genes");

        } else {
            hlaGenes = null;
        }


        

        // Split the null GWAS matrix into the parts used in pathway enrichment.
        final int nrSampleToUseForCorrelation = options.getPermutationGeneCorrelations();
        final int nrSamplesToUseForNullBetas = options.getPermutationPathwayEnrichment();
        final int nrSamplesToUseForFDR = options.getPermutationFDR();

        final Set<String> nullGwasRuns = genePvaluesNullGwas.getHashCols().keySet();
        if (nullGwasRuns.size() < (nrSampleToUseForCorrelation + nrSamplesToUseForNullBetas)) {
            throw new Exception("Not enough null gwas runs: " + nullGwasRuns.size() + " < " + nrSampleToUseForCorrelation + " + " + nrSamplesToUseForNullBetas + " + " + nrSamplesToUseForFDR);
        }

        Iterator<String> nullGwasRunIterator = nullGwasRuns.iterator();

        final LinkedHashSet<String> sampleToUseForCorrelation = new LinkedHashSet<>(nrSampleToUseForCorrelation);
        for (int i = 0; i < nrSampleToUseForCorrelation; ++i) {
            sampleToUseForCorrelation.add(nullGwasRunIterator.next());
        }

        // Null beta and FDR samples are combined in the same matrix and split later in PathwayEnrichments as
        // beta calculation is independent of other columns
        final LinkedHashSet<String> samplesToUseForNullBetas = new LinkedHashSet<>(nrSamplesToUseForNullBetas + nrSamplesToUseForFDR);
        for (int i = 0; i < nrSamplesToUseForNullBetas + nrSamplesToUseForFDR; ++i) {
            samplesToUseForNullBetas.add(nullGwasRunIterator.next());
        }

        final DoubleMatrixDataset<String, String> geneZscoresNullGwasCorrelation = genePvaluesNullGwas.viewColSelection(sampleToUseForCorrelation);
        final DoubleMatrixDataset<String, String> geneZscoresNullGwasNullBetas = genePvaluesNullGwas.viewColSelection(samplesToUseForNullBetas);

        final DoubleMatrixDataset<String, String> geneMaxSnpZscoreNullGwasCorrelation = geneMaxSnpZscoreNullGwas.viewColSelection(sampleToUseForCorrelation);
        final DoubleMatrixDataset<String, String> geneMaxSnpZscoreNullGwasBetas = geneMaxSnpZscoreNullGwas.viewColSelection(samplesToUseForNullBetas);

        ArrayList<PathwayEnrichments> pathwayEnrichments = new ArrayList<>(pathwayDatabases.size());
        for (PathwayDatabase pathwayDatabase : pathwayDatabases) {
            pathwayEnrichments.add(new PathwayEnrichments(
                    pathwayDatabase,
                    selectedGenes,
                    genes,
                    options.isForceNormalPathwayPvalues(),
                    options.isForceNormalGenePvalues(),
                    genePvalues,
                    geneZscoresNullGwasCorrelation,
                    geneZscoresNullGwasNullBetas,
                    options.getOutputBasePath(),
                    hlaGenes,
                    options.isIgnoreGeneCorrelations(),
                    options.getGenePruningR(),
                    options.getDebugFolder(),
                    options.getIntermediateFolder(),
                    options.isQuantileNormalizePermutations(),
                    options.isRegressGeneLengths(),
                    geneMaxSnpZscore,
                    geneMaxSnpZscoreNullGwasCorrelation,
                    geneMaxSnpZscoreNullGwasBetas,
                    options.getGeneCorrelationWindow(),
                    nrSamplesToUseForFDR
            ));
        }

        LOGGER.info("Completed enrichment analysis for " + pathwayDatabases.size() + " pathway databases");
        return new DownstreamerStep2Results(pathwayEnrichments, genePvalues);
    }

    /**
     * Assign genes to GWAS loci based on a given window, will be used for cis-priotiziation of GWAS genes
     *
     * @param options
     * @return DownstreamerStep3Results
     * @throws Exception
     */
    public static DownstreamerStep3Results step3(DownstreamerOptions options) throws Exception {
        // TODO: paramterize
        int window = 500000;

        // GWAS pvalues
        DoubleMatrixDataset<String, String> gwasSnpZscores = DoubleMatrixDataset.loadDoubleBinaryData(options.getGwasZscoreMatrixPath());

        // Gene info
        Map<String, List<Gene>> genes = IoUtils.readGenesAsChrMap(options.getGeneInfoFile());
        LOGGER.info("Loaded " + genes.size() + " genes");

        // Independent variants, defined in step1
        Map<String, Set<String>> independentVariants = IoUtils.readIndependentVariants(options.getGwasTopHitsFile());
        Map<String, Set<String>> alternativeIndependentVariants = IoUtils.readAlternativeIndependentVariants(options.getAlternativeTopHitFiles());

        // Optional alternative top hits, will override internal depict prio
        Map<String, File> alternativeTopHitFiles = options.getAlternativeTopHitFiles();

        // Determine the set of possible variants where info needs to be retrieved
        Set<String> possibleIndexVariants = new HashSet<>();
        for (String trait : gwasSnpZscores.getColObjects()) {
            possibleIndexVariants.addAll(independentVariants.get(trait));

            if (alternativeTopHitFiles.containsKey(trait)) {
                possibleIndexVariants.addAll(alternativeIndependentVariants.get(trait));
            }
        }

        // Genotype data (for positions, not doing ld clumping here)
        RandomAccessGenotypeData referenceGenotypeData = IoUtils.readReferenceGenotypeDataMatchingGwasSnps(options, possibleIndexVariants);
        Map<String, GeneticVariant> variantMap = referenceGenotypeData.getVariantIdMap();
        LOGGER.info("Done loading genotype data");

        // Output store
        DownstreamerStep3Results output = new DownstreamerStep3Results();

        for (String trait : gwasSnpZscores.getColObjects()) {

            Map<String, SummaryStatisticRecord> records = new HashMap<>();
            List<SummaryStatisticRecord> recordCache = new ArrayList<>();

            // Select index variants and overlap with reference genotypes
            Set<String> curIndexVariants = independentVariants.get(trait);
            curIndexVariants.retainAll(variantMap.keySet());

            // Use either internal Downstreamer topsnps or the ones from the  provided file
            if (alternativeTopHitFiles.containsKey(trait)) {
                LOGGER.info("Using alternative top hits for: " + trait);
                recordCache = IoUtils.readAlternativeIndependentVariantsAsRecords(alternativeTopHitFiles.get(trait));
            } else {
                for (String variant : curIndexVariants) {
                    double pvalue = ZScores.zToP(gwasSnpZscores.getElement(variant, trait));
                    recordCache.add(new SummaryStatisticRecord(variantMap.get(variant), pvalue));
                }
            }
            LOGGER.info("Parsed " + recordCache.size() + " index variants (independent SNPs)");

            // Convert cache to SNP id Map
            for (SummaryStatisticRecord outRec : recordCache) {
                records.put(outRec.getPrimaryVariantId(), outRec);
            }

            // Generate loci based on predefined index variants
            List<Locus> loci = LocusUtils.makeLociWithGivenIndexSnps(records, curIndexVariants, window);
            LOGGER.info("Made " + loci.size() + " loci (independent SNPs) using " + window + " bp window");

            for (Locus curLocus : loci) {
                // Annotate the independent top variants
                if (alternativeTopHitFiles.containsKey(trait)) {
                    curLocus.addIndepVariants(alternativeIndependentVariants.get(trait));
                } else {
                    curLocus.addIndepVariants(independentVariants.get(trait));
                }

                // Add the genes to the loci objects
                // TODO: not very efficient, but should be fine, the existing implementations were not really easily implementable. Alternatively could look into using R-trees if this proves to slow
                if (genes.containsKey(curLocus.getSequenceName())) {
                    for (Gene curGene : genes.get(curLocus.getSequenceName())) {
                        if (curGene.isOverlapping(curLocus)) {
                            curLocus.addGene(curGene);
                        }
                    }
                } else {
                    LOGGER.warn("Omitting " + curLocus.getSequenceName() + ":" + curLocus.getStart() + "-" + curLocus.getEnd() + " as the chromosome could not be found");
                }
            }

            output.addLoci(trait, loci);
        }


        return output;
    }


    private static double[] generateRandomChi2(long numberOfPermutations, int numberOfVariantPerGeneToExpect) {

        final double[] randomChi2;
        if ((numberOfPermutations * numberOfVariantPerGeneToExpect) > Integer.MAX_VALUE - 10) {
            randomChi2 = new double[Integer.MAX_VALUE - 10];
        } else {
            randomChi2 = new double[numberOfVariantPerGeneToExpect * (int) numberOfPermutations];
        }

        final int randomChi2Size = randomChi2.length;
        final int nrThreads = DownstreamerOptions.getNumberOfThreadsToUse();
        final int permPerThread = randomChi2Size / nrThreads;
        final int leftoverPerm = randomChi2Size % nrThreads;

        IntStream.range(0, nrThreads).parallel().forEach(task -> {

            final ThreadLocalRandom rnd = ThreadLocalRandom.current();
            double z;
            for (int p = 0; p < permPerThread; ++p) {
                z = rnd.nextGaussian();
                randomChi2[(task * permPerThread) + p] = z * z;
            }

        });

        if (leftoverPerm > 0) {
            final ThreadLocalRandom rnd = ThreadLocalRandom.current();
            double z;
            for (int p = 0; p < leftoverPerm; ++p) {
                z = rnd.nextGaussian();
                randomChi2[(nrThreads * permPerThread) + p] = z * z;
            }
        }

        return randomChi2;

    }
}
