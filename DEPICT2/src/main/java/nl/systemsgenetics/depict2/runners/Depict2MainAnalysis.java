package nl.systemsgenetics.depict2.runners;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import nl.systemsgenetics.depict2.*;
import nl.systemsgenetics.depict2.gene.Gene;
import nl.systemsgenetics.depict2.gene.GenePvalueCalculator;
import nl.systemsgenetics.depict2.io.IoUtils;
import nl.systemsgenetics.depict2.pathway.PathwayDatabase;
import nl.systemsgenetics.depict2.pathway.PathwayEnrichments;
import nl.systemsgenetics.depict2.summarystatistic.Locus;
import nl.systemsgenetics.depict2.summarystatistic.LocusUtils;
import nl.systemsgenetics.depict2.summarystatistic.SummaryStatisticRecord;
import nl.systemsgenetics.depict2.summarystatistic.filters.PvalueFilterSmaller;
import org.apache.log4j.Logger;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.ZScores;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ThreadLocalRandom;
import java.util.stream.IntStream;

/**
 * Runners for the main Depict 2 analysis.
 */
public class Depict2MainAnalysis {

    private static final Logger LOGGER = Logger.getLogger(Depict2MainAnalysis.class);

    /**
     * Calculate the real gene p-values from GWAS z-scores as well as those for x ammount of random GWAS phenotypes.
     *
     * @param options
     * @return Depict2Step1Results
     * @throws IOException
     * @throws Exception
     */
    public static Depict2Step1Results step1(Depict2Options options) throws IOException, Exception {

        //Test here to prevent chrash after running for a while
        if (!new File(options.getGwasZscoreMatrixPath() + ".dat").exists()) {
            throw new FileNotFoundException("GWAS matrix does not exist at: " + options.getGwasZscoreMatrixPath() + ".dat");
        }

        RandomAccessGenotypeData referenceGenotypeData = readReferenceGenotypeDataMatchingGwasSnps(options);
        LOGGER.info("Done loading genotype data");

        List<Gene> genes = IoUtils.readGenes(options.getGeneInfoFile());
        LOGGER.info("Loaded " + genes.size() + " genes");

        double[] randomChi2 = generateRandomChi2(options.getNumberOfPermutationsRescue(), 500);

        LOGGER.info("Prepared reference null distribution with " + Depict2.LARGE_INT_FORMAT.format(randomChi2.length) + " values");

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

        return new Depict2Step1Results(genePvalues, genePvaluesNullGwas, geneVariantCount, geneMaxSnpZscore, geneMaxSnpZscoreNullGwas);
    }

    /**
     * Calculate the association between gene pvalues and pathway databases using a GLS model that accounts
     * for the correlation between gene-pvalues.
     *
     * @param options
     * @param step1Res
     * @return Depict2Step2Results
     * @throws IOException
     * @throws Exception
     */
    public static Depict2Step2Results step2(Depict2Options options, Depict2Step1Results step1Res) throws IOException, Exception {

        options.getIntermediateFolder().mkdir();

        if (options.getMode() == Depict2Mode.STEP2) {
            LOGGER.info("Continuing previous analysis by loading gene p-values");
            step1Res = Depict2Step1Results.loadFromDisk(options.getRun1BasePath());
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

        LOGGER.info("Number of genes with atleast one variant in specified window: " + Depict2.LARGE_INT_FORMAT.format(selectedGenes.size()));
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


        final List<PathwayDatabase> pathwayDatabases = options.getPathwayDatabases();

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
        return new Depict2Step2Results(pathwayEnrichments, genePvalues);
    }

    /**
     * Assign genes to GWAS loci based on a given window, will be used for cis-priotiziation of GWAS genes
     *
     * @param options
     * @return Depict2Step3Results
     * @throws Exception
     */
    public static Depict2Step3Results step3(Depict2Options options) throws Exception {
        // TODO: paramterize
        double upperPvalThresh = 5e-8;
        double lowerPvalThresh = upperPvalThresh;
        int window = 1000000;

        // GWAS pvalues
        DoubleMatrixDataset<String, String> gwasSnpZscores = DoubleMatrixDataset.loadDoubleBinaryData(options.getGwasZscoreMatrixPath());

        // Filter variants on pvalue
        Set<String> significantVariants = new HashSet<>();
        for (String trait : gwasSnpZscores.getColObjects()) {
            for (String variant : gwasSnpZscores.getColObjects()) {
                double pvalue = ZScores.zToP(gwasSnpZscores.getElement(variant, trait));
                if (pvalue < upperPvalThresh) {
                    significantVariants.add(variant);
                }
            }
        }
        gwasSnpZscores = gwasSnpZscores.viewRowSelection(significantVariants);
        LOGGER.info("Done loading GWAS data, " + significantVariants + " variants at p < " + upperPvalThresh);

        // Genotype data (for positions, not doing ld clumping here)
        RandomAccessGenotypeData referenceGenotypeData = readReferenceGenotypeDataMatchingGwasSnps(options, significantVariants);
        Map<String, GeneticVariant> variantMap = referenceGenotypeData.getVariantIdMap();
        LOGGER.info("Done loading genotype data");

        // Gene info
        Map<String, List<Gene>> genes = IoUtils.readGenesAsChrMap(options.getGeneInfoFile());
        LOGGER.info("Loaded " + genes.size() + " genes");

        // Independent variants, defined in step1
        Map<String, Set<String>> independentVariants = IoUtils.readIndependentVariants(options.getOutputBasePath() + "_independentTopVariants.txt");

        // Output store
        Depict2Step3Results output = new Depict2Step3Results();

        for (String trait : gwasSnpZscores.getColObjects()) {
            Map<String, SummaryStatisticRecord> records = new HashMap<>();
            for (String variant : gwasSnpZscores.getColObjects()) {
                double pvalue = ZScores.zToP(gwasSnpZscores.getElement(variant, trait));
                if (pvalue < upperPvalThresh) {
                    records.put(variant, new SummaryStatisticRecord(variantMap.get(variant), pvalue));
                }
            }

            // Generate loci based on GWAS topsnps
            List<Locus> loci = LocusUtils.makeLoci(records, new PvalueFilterSmaller(upperPvalThresh), new PvalueFilterSmaller(lowerPvalThresh), window);
            LOGGER.info("Detected " + loci.size() + " loci using " + window + " bp window");

            for (Locus curLocus : loci) {
                // Annotate the independent top variants
                curLocus.addIndepVariants(independentVariants.get(trait));

                // Add the genes to the loci objects
                // TODO: not very efficient, but should be fine, the existing implementations were not really easily implementable. Alternatively could look into using R-trees if this proves to slow
                for (Gene curGene : genes.get(curLocus.getSequenceName())) {
                    if (curGene.isOverlapping(curLocus)) {
                        curLocus.addGene(curGene);
                    }
                }
            }

            output.addLoci(trait, loci);
        }

        return output;
    }


    /**
     * Load genotype data matching GWAS matrix and MAF filter.
     *
     * @param options Depict options object
     * @return RandomAccesGenotypeData for all SNPs in GWAS matrix and MAF
     * @throws IOException
     */
    private static RandomAccessGenotypeData readReferenceGenotypeDataMatchingGwasSnps(Depict2Options options) throws IOException {
        return readReferenceGenotypeDataMatchingGwasSnps(options, null);
    }


    /**
     * Load genotype data matching GWAS matrix and MAF filter.
     *
     * @param options Depict options object
     * @return RandomAccesGenotypeData for all SNPs in GWAS matrix and MAF
     * @throws IOException
     */
    private static RandomAccessGenotypeData readReferenceGenotypeDataMatchingGwasSnps(Depict2Options options, Set<String> variantSubset) throws IOException {


        final List<String> variantsInZscoreMatrix;
        if (variantSubset == null) {
            variantsInZscoreMatrix = IoUtils.readMatrixAnnotations(new File(options.getGwasZscoreMatrixPath() + ".rows.txt"));
        } else {
            variantsInZscoreMatrix = new ArrayList<>(variantSubset);
        }

        final List<String> phenotypesInZscoreMatrix = IoUtils.readMatrixAnnotations(new File(options.getGwasZscoreMatrixPath() + ".cols.txt"));

        LOGGER.info("Number of phenotypes in GWAS matrix: " + Depict2.LARGE_INT_FORMAT.format(phenotypesInZscoreMatrix.size()));
        LOGGER.info("Number of variants in GWAS matrix: " + Depict2.LARGE_INT_FORMAT.format(variantsInZscoreMatrix.size()));

        if (options.getVariantFilterFile() != null) {
            HashSet<String> variantsToInclude = IoUtils.readVariantFilterFile(options.getVariantFilterFile());
            Iterator<String> variantsInZscoreMatrixIt = variantsInZscoreMatrix.iterator();
            while (variantsInZscoreMatrixIt.hasNext()) {
                String variant = variantsInZscoreMatrixIt.next();
                if (!variantsToInclude.contains(variant)) {
                    variantsInZscoreMatrixIt.remove();
                }
            }
            LOGGER.info("Number of variants after filtering on selected variants: " + Depict2.LARGE_INT_FORMAT.format(variantsInZscoreMatrix.size()));
        }

        return IoUtils.loadGenotypes(options, variantsInZscoreMatrix);
    }


    private static double[] generateRandomChi2(long numberOfPermutations, int numberOfVariantPerGeneToExpect) {

        final double[] randomChi2;
        if ((numberOfPermutations * numberOfVariantPerGeneToExpect) > Integer.MAX_VALUE - 10) {
            randomChi2 = new double[Integer.MAX_VALUE - 10];
        } else {
            randomChi2 = new double[numberOfVariantPerGeneToExpect * (int) numberOfPermutations];
        }

        final int randomChi2Size = randomChi2.length;
        final int nrThreads = Depict2Options.getNumberOfThreadsToUse();
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
