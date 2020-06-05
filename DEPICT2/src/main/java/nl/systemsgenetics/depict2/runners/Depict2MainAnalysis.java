package nl.systemsgenetics.depict2.runners;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import nl.systemsgenetics.depict2.Depict2;
import nl.systemsgenetics.depict2.Depict2Mode;
import nl.systemsgenetics.depict2.Depict2Options;
import nl.systemsgenetics.depict2.gene.Gene;
import nl.systemsgenetics.depict2.gene.GenePvalueCalculator;
import nl.systemsgenetics.depict2.io.ExcelWriter;
import nl.systemsgenetics.depict2.io.IoUtils;
import nl.systemsgenetics.depict2.pathway.PathwayDatabase;
import nl.systemsgenetics.depict2.pathway.PathwayEnrichments;
import org.apache.log4j.Logger;
import org.molgenis.genotype.RandomAccessGenotypeData;
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
 *
 */
public class Depict2MainAnalysis {
    
    private static final Logger LOGGER = Logger.getLogger(Depict2MainAnalysis.class);

    public static void run(Depict2Options options) throws IOException, Exception {

        //Test here to prevent chrash after running for a while
        if (!new File(options.getGwasZscoreMatrixPath() + ".dat").exists()) {
            throw new FileNotFoundException("GWAS matrix does not exist at: " + options.getGwasZscoreMatrixPath() + ".dat");
        }

        final List<String> variantsInZscoreMatrix = IoUtils.readMatrixAnnotations(new File(options.getGwasZscoreMatrixPath() + ".rows.txt"));
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

        RandomAccessGenotypeData referenceGenotypeData = IoUtils.loadGenotypes(options, variantsInZscoreMatrix);

        LOGGER.info("Done loading genotype data");

        List<Gene> genes = IoUtils.readGenes(options.getGeneInfoFile());

        LOGGER.info("Loaded " + genes.size() + " genes");

        double[] randomChi2 = generateRandomChi2(options.getNumberOfPermutationsRescue(), 500);

        LOGGER.info("Prepared reference null distribution with " + Depict2.LARGE_INT_FORMAT.format(randomChi2.length) + " values");

        File usedVariantsPerGeneFile = options.isSaveUsedVariantsPerGene() ? new File(options.getOutputBasePath() + "_usedVariantsPerGene.txt") : null;

        GenePvalueCalculator gpc = new GenePvalueCalculator(options.getGwasZscoreMatrixPath(), referenceGenotypeData, genes, options.getWindowExtend(), options.getMaxRBetweenVariants(), options.getNumberOfPermutations(), options.getNumberOfPermutationsRescue(), options.getOutputBasePath(), randomChi2, options.correctForLambdaInflation(), options.getPermutationGeneCorrelations(), options.getPermutationPathwayEnrichment(), options.getDebugFolder(), options.getVariantGeneLinkingFile(), usedVariantsPerGeneFile);

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
        LOGGER.info("Gene p-values saved. If needed the analysis can be resummed from this point using --mode RUN2 and exactly the same output path and genes file");

        if (options.getPathwayDatabases().isEmpty()) {
            LOGGER.info("The analysis will now stop since no pathway databases are provided. Use --mode RUN2 and exactly the same output path and genes file to continue");
        } else {
            run2(options, genePvalues, genePvaluesNullGwas, genes, geneVariantCount, geneMaxSnpZscore, geneMaxSnpZscoreNullGwas);
        }
    }

    /**
     * @param options
     * @param genePvalues
     * @param genePvaluesNullGwas
     * @param genes
     * @throws IOException
     * @throws Exception
     */
    public static void run2(Depict2Options options, DoubleMatrixDataset<String, String> genePvalues, DoubleMatrixDataset<String, String> genePvaluesNullGwas, List<Gene> genes, DoubleMatrixDataset<String, String> geneVariantCount, DoubleMatrixDataset<String, String> geneMaxSnpZscore, DoubleMatrixDataset<String, String> geneMaxSnpZscoreNullGwas) throws IOException, Exception {

        options.getIntermediateFolder().mkdir();

        if (options.getMode() == Depict2Mode.STEP2) {
            LOGGER.info("Continuing previous analysis by loading gene p-values");
            if (new File(options.getRun1BasePath()+ "_genePvalues.dat").exists()) {
                genePvalues = DoubleMatrixDataset.loadDoubleBinaryData(options.getRun1BasePath() + "_genePvalues");
                genePvaluesNullGwas = DoubleMatrixDataset.loadDoubleBinaryData(options.getRun1BasePath() + "_genePvaluesNullGwas");

                // Always load to avoid nullpointers
                geneMaxSnpZscore = DoubleMatrixDataset.loadDoubleBinaryData(options.getRun1BasePath() + "_geneMaxSnpScores");
                geneMaxSnpZscoreNullGwas = DoubleMatrixDataset.loadDoubleBinaryData(options.getRun1BasePath() + "_geneMaxSnpZscoresNullGwas");

            } else {
                LOGGER.fatal("Could not find gene pvalues at: " + options.getRun1BasePath() + "_genePvalues.dat");
                LOGGER.fatal("First use --mode RUN to calculate gene p-values");
                return;
            }
            geneVariantCount = DoubleMatrixDataset.loadDoubleTextData(options.getRun1BasePath() + "_geneVariantCount.txt", '\t');
            LOGGER.info("Gene p-values loaded");
            genes = IoUtils.readGenes(options.getGeneInfoFile());
            LOGGER.info("Loaded " + genes.size() + " genes");
        }

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
            for (Gene gene : genes) {
                if (gene.getChr().equals("6") && ((gene.getStart() > 20000000 && gene.getStart() < 40000000) || (gene.getStop() > 20000000 && gene.getStop() < 40000000))) {
                    hlaGenes.add(gene.getGene());
                }
            }
            LOGGER.info("Excluding " + hlaGenes.size() + " genes");

        } else {
            hlaGenes = null;
        }

        final List<PathwayDatabase> pathwayDatabases = options.getPathwayDatabases();
        final int nrSampleToUseForCorrelation = options.getPermutationGeneCorrelations();
        final int nrSamplesToUseForNullBetas = options.getPermutationPathwayEnrichment();

        final Set<String> nullGwasRuns = genePvaluesNullGwas.getHashCols().keySet();
        if (nullGwasRuns.size() < (nrSampleToUseForCorrelation + nrSamplesToUseForNullBetas)) {
            throw new Exception("Not enough null gwas runs: " + nullGwasRuns.size() + " < " + nrSampleToUseForCorrelation + " + " + nrSamplesToUseForNullBetas);
        }

        Iterator<String> nullGwasRunIterator = nullGwasRuns.iterator();

        final LinkedHashSet<String> sampleToUseForCorrelation = new LinkedHashSet<>(nrSampleToUseForCorrelation);
        for (int i = 0; i < nrSampleToUseForCorrelation; ++i) {
            sampleToUseForCorrelation.add(nullGwasRunIterator.next());
        }

        final LinkedHashSet<String> samplesToUseForNullBetas = new LinkedHashSet<>(nrSamplesToUseForNullBetas);
        for (int i = 0; i < nrSamplesToUseForNullBetas; ++i) {
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
                    geneMaxSnpZscoreNullGwasBetas
            ));
        }

        if (options.isSaveOuputAsExcelFiles()) {
            ExcelWriter.saveEnrichmentsToExcel(pathwayEnrichments, options.getOutputBasePath(), genePvalues.getColObjects(), hlaGenes != null,options);
        } else {
            for (PathwayEnrichments pathwayEnrichment : pathwayEnrichments) {
                //this will make sure z-scores are saved even if make excel is off
                pathwayEnrichment.getEnrichmentZscores();
            }
        }

        LOGGER.info("Completed enrichment analysis for " + pathwayDatabases.size() + " pathway databases");

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
