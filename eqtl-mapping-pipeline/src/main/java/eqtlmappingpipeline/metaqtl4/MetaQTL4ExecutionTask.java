/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl4;

import gnu.trove.map.hash.TObjectIntHashMap;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletionService;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.util.Primitives;

/**
 *
 * @author Harm-Jan
 */
class MetaQTL4ExecutionTask implements Callable<Boolean> {

    private final long[] randomizationSeeds;
    private final ArrayList<MetaQTL4MetaTrait> availableTraits;
    private final TObjectIntHashMap<MetaQTL4MetaTrait> availableTraitsHash;
    private final MetaQTL4Dataset[] datasets;
    private final GeneticVariant[][] geneticVariantIndex;
    private final MetaQTL4Settings m_settings;
    private final MetaQTL4TraitAnnotation traitAnnotation;
    private final Integer[][] traitIndex;
    private final Set<MetaQTL4MetaTrait> traitsToInclude;
    private final Set<String> variantsToInclude;
    private final int bufferSize;
    private final int start;
    private final int stop;
    private final CompletionService resultPool;

    public MetaQTL4ExecutionTask(int nrThreads, long[] randomizationSeeds, ArrayList<MetaQTL4MetaTrait> availableTraits, TObjectIntHashMap<MetaQTL4MetaTrait> availableTraitsHash, MetaQTL4Dataset[] datasets, GeneticVariant[][] geneticVariantIndex, MetaQTL4Settings m_settings, MetaQTL4TraitAnnotation traitAnnotation, Integer[][] traitIndex, Set<MetaQTL4MetaTrait> traitsToInclude, Set<String> variantsToInclude, int start, int stop, int bufferSize, CompletionService resultPool) {

        this.randomizationSeeds = randomizationSeeds;
        this.availableTraits = availableTraits;
        this.availableTraitsHash = availableTraitsHash;
        this.datasets = datasets;
        this.geneticVariantIndex = geneticVariantIndex;
        this.m_settings = m_settings;
        this.traitAnnotation = traitAnnotation;
        this.traitIndex = traitIndex;
        this.traitsToInclude = traitsToInclude;
        this.variantsToInclude = variantsToInclude;
        this.start = start;
        this.stop = stop;
        this.bufferSize = bufferSize;
        this.resultPool = resultPool;
    }

    @Override
    public Boolean call() throws Exception {
        try {
            // normal procedure
            int nrPermutations = randomizationSeeds.length;
            int nrDatasets = geneticVariantIndex.length;
            int nrVariants = geneticVariantIndex[nrDatasets - 1].length;

            // iterate the variants
            int numberOfDatasetsPassingQC = 0;

            double mafthreshold = m_settings.getSnpQCMAFThreshold();
            double hwepthreshold = m_settings.getSnpQCHWEThreshold();
            double callratethreshold = m_settings.getSnpQCHWEThreshold();

            int[] maxNrSamplesPerDataset = new int[datasets.length];
            boolean[][] includeTraitSample = new boolean[nrDatasets][0];
            for (int d = 0; d < maxNrSamplesPerDataset.length; d++) {
                int nrSamples = datasets[d].getGenotypeToTraitCouplingInt().length;
                maxNrSamplesPerDataset[d] = nrSamples;
                includeTraitSample[d] = new boolean[nrSamples];
            }

            ProgressBar pb = new ProgressBar(nrVariants, "Running calculations...");

            double[][][] genotypes = new double[bufferSize][nrDatasets][0];
            double[][] genotypeVariances = new double[bufferSize][nrDatasets];
            double[][] sampleSizes = new double[bufferSize][nrDatasets];
            boolean[][] flipEffects = new boolean[bufferSize][nrDatasets];
            GeneticVariant[][] variantBuffer = new GeneticVariant[bufferSize][nrDatasets];
            int lastLoadedVariant = start;
            int nrToLoad = bufferSize;
            while (lastLoadedVariant < stop) {

                if (lastLoadedVariant + bufferSize >= stop) {
                    nrToLoad = stop - lastLoadedVariant;
                }
                for (int datasetId = 0; datasetId < nrDatasets; datasetId++) {
                    for (int bufferPosition = 0; bufferPosition < nrToLoad; bufferPosition++) {
                        int indexPosition = lastLoadedVariant + bufferPosition;

                        GeneticVariant variant = geneticVariantIndex[datasetId][indexPosition];

                        MetaQTL4Dataset dataset = datasets[datasetId];
                        int[] gte = dataset.getGenotypeToTraitCouplingInt();

                        if (variant != null) {
                            float[] genotypesTMP = dataset.getSampleDosages(variant);
                            double maf = variant.getMinorAlleleFrequency();
                            double hwep = variant.getHwePvalue();
                            double callrate = variant.getCallRate();

                            // TODO: best done on loading dataset with filter.
                            if (maf > mafthreshold && hwep > hwepthreshold && callrate > callratethreshold) {
                                // TODO: remove missing genotypes and rescale the genotype data
                                Pair<double[], Double> genotypedata = correctGenotypesForMissingValuesAndNormalize(gte, variant, genotypesTMP, includeTraitSample[datasetId]);
                                genotypes[bufferPosition][datasetId] = genotypedata.getLeft();

                                genotypeVariances[bufferPosition][datasetId] = genotypedata.getRight();
                                sampleSizes[bufferPosition][datasetId] = genotypes[datasetId].length;
                                variantBuffer[bufferPosition][datasetId] = variant;
                            } else {
                                // TODO: output to the log?
                                genotypes[bufferPosition][datasetId] = null;
                                genotypeVariances[bufferPosition][datasetId] = 0;
                                variantBuffer[bufferPosition][datasetId] = null;
                            }
                        }
                    }
                }

                // SNP buffer loaded..
                // determine list of probes for each variant
                ArrayList<HashSet<Integer>> traitsToTest = null; //
                for (int bufferPosition = 0; bufferPosition < nrToLoad; bufferPosition++) {
                    int startPosition = -1;
                    String sequence = null;
                    for (int dataset = 0; dataset < datasets.length; dataset++) {
                        GeneticVariant variant = variantBuffer[bufferPosition][dataset];
                        if (sequence == null) {
                            startPosition = variant.getStartPos();
                            sequence = variant.getSequenceName();
                        }
                    }
                    if (m_settings.isCisAnalysis() ^ m_settings.isTransAnalysis()) {
                        // index once, use many
                        if (traitsToTest == null) {
                            traitsToTest = new ArrayList<HashSet<Integer>>();
                        }
                        HashSet<Integer> traitsToTestTMP = new HashSet<Integer>();
                        // possible faster with hashmap lookup
                        for (int i = 0; i < availableTraits.size(); i++) {
                            MetaQTL4MetaTrait trait = availableTraits.get(i); // a treemap may be more efficient here
                            boolean sameChr = trait.getChr().equals(sequence);
                            if (sameChr) {
                                int distance = Math.abs(trait.getChrMidpoint() - startPosition); // precalculate this, if possible
                                if (m_settings.isCisAnalysis() && distance < m_settings.getCiseQTLAnalysMaxSNPProbeMidPointDistance()) {
                                    traitsToTestTMP.add(i);
                                }
                                if (m_settings.isTransAnalysis() && distance < m_settings.getCiseQTLAnalysMaxSNPProbeMidPointDistance()) {
                                    traitsToTestTMP.add(i);
                                }
                            }
                        }
                        traitsToTest.add(traitsToTestTMP);
                    }
                }

                // run the analysis
                for (int permutation = -1; permutation < nrPermutations; permutation++) {
                    // iterate the variants per dataset..
                    for (int bufferPosition = 0; bufferPosition < nrToLoad; bufferPosition++) {
                        double[][] datasetGenotypes = genotypes[bufferPosition];

                        // permute the genotypes, if required.
                        if (permutation >= 0) {
                            long seed = randomizationSeeds[permutation];
                            Random rand = new Random(seed);
                            // reorder genotypes on the basis of the permuted sample links
                            for (int d = 0; d < datasets.length; d++) {
                                if (datasetGenotypes[d] != null) {
                                    datasetGenotypes[d] = permuteGenotypes(datasetGenotypes[d], rand);
                                }
                            }
                        }

                        if (traitsToTest == null) {
                            // cistrans analysis
                            for (int trait = 0; trait < availableTraits.size(); trait++) {
                                //int bin = test(sampleSizes, genotypes, genotypeVariances, trait, includeTraitSample);

                            }
                        } else {
                            HashSet<Integer> traitsToTestForSNP = traitsToTest.get(bufferPosition);
                            if (m_settings.isCisAnalysis()) {

                                for (Integer trait : traitsToTestForSNP) {
//                                MetaQTL4MetaTrait traitObj = availableTraits.get(trait);
//                                System.out.println(Strings.concat(traitObj.getPlatformIds(), Strings.comma));
                                    //int bin = test(sampleSizes, genotypes, genotypeVariances, trait, includeTraitSample);

                                }
                            }
                            if (m_settings.isTransAnalysis()) {
                                for (int i = 0; i < availableTraits.size(); i++) {
                                    if (!traitsToTestForSNP.contains(i)) {

                                    }
                                }

                            }

                        }

                    }

                }

                lastLoadedVariant += nrToLoad;
            }

            pb.close();

        } catch (Exception e) {
            e.printStackTrace();
        }
        return true;
    }

    private void test(int[] sampleSizes, double[][] genotypes, double[] genotypeVariances, Integer traitId, boolean[][] includeTraitSample) {
        int nrDatasets = datasets.length;
        double[] zscores = new double[nrDatasets];
        for (int datasetId = 0; datasetId < nrDatasets; datasetId++) {

            if (genotypes[datasetId] == null) {
                zscores[datasetId] = Double.NaN;
            } else {
                double[] x = genotypes[datasetId];
                if (x == null) {
                    System.err.println("ERROR: genotype is null");
                    zscores[datasetId] = Double.NaN;
                } else {
                    double varianceX = genotypeVariances[datasetId];
                    boolean[] includeDatasetTraitSample = includeTraitSample[datasetId];

                    double meanY;
                    double varianceY;

                    // re-normalize the trait data when the genotypes had missing values
                    Integer datasetTraitId = traitIndex[datasetId][traitId];
                    double[] y = datasets[datasetId].getTraitData(datasetTraitId);
//                    if (sampleSizes[datasetId] != y.length) {
                    double[] newY = new double[x.length];
                    int itr = 0;
                    double sum = 0;

                    // recalculate mean and variance
                    for (int s = 0; s < y.length; s++) {
//                        if (includeDatasetTraitSample[s]) {
                        newY[itr] = y[s];
                        sum += newY[itr];
                        itr++;
//                        }
                    }

                    y = newY;
                    meanY = sum / itr;
                    double varsum = 0;
                    for (int i = 0; i < y.length; i++) {
                        y[i] -= meanY;
                        varsum += y[i] * y[i];
                    }
                    varianceY = varsum / (y.length - 1);
//                    } else {
                    double varianceYOrig = datasets[datasetId].getTraitVariance(traitId);
//                    }

                    if (varianceY == 0) {
                        // trait has no variance
                        zscores[datasetId] = Double.NaN;
                    } else {
                        //Calculate correlation coefficient:
                        double correlation = Correlation.correlateMeanCenteredData(x, y, varianceX, varianceY);
                        double correlation2 = JSci.maths.ArrayMath.correlation(x, y);
                        System.out.println(correlation + "\t" + correlation2 + "\t" + x.length + "\t" + varianceX + "\t" + varianceYOrig + "\t" + varianceY + "\t" + Descriptives.variance(x) + "\t" + Descriptives.variance(y));
                        if (correlation >= -1 && correlation <= 1) {
                            zscores[datasetId] = Correlation.convertCorrelationToZScore(x.length, correlation);
                        } else {
                            System.err.println("Error! correlation invalid: " + correlation);
                            System.exit(-1);
                        }
                    }
                }
            }
        }

        double metaZ = ZScores.getWeightedZ(zscores, sampleSizes);
        double p = ZScores.zToP(metaZ);

    }

    private Pair<double[], Double> correctGenotypesForMissingValuesAndNormalize(int[] coupledTraitInt, GeneticVariant variant, float[] genotypesTMP, boolean[] includeTraitSample) {
        int xLen = genotypesTMP.length;
        int nrMissing = 0;
        double[] snpmeancorrectedgenotypes;
        double varianceX;
        double meanX;
        if (variant.getCallRate() < 1d) {
            double sum = 0;
            for (int i = 0; i < xLen; i++) {
                if (genotypesTMP[i] < 0) {
                    nrMissing++;
                    int coupledTrait = coupledTraitInt[i]; // not really sure if this is still required..
                    includeTraitSample[coupledTrait] = false;
                } else {
                    sum += genotypesTMP[i];
                }
            }
            int newXLen = xLen - nrMissing;
            meanX = sum / newXLen;
            snpmeancorrectedgenotypes = new double[newXLen];

            int itr = 0;
            for (int i = 0; i < newXLen; i++) {
                if (genotypesTMP[i] >= 0) {
                    snpmeancorrectedgenotypes[itr] = genotypesTMP[i] - meanX;
                    itr++;
                }
            }
        } else {
            double sum = 0;
            snpmeancorrectedgenotypes = new double[xLen];
            for (int i = 0; i < xLen; i++) {
                double genotype = genotypesTMP[i];
                sum += genotype;
                snpmeancorrectedgenotypes[i] = genotype;
            }
            //  meanX = sum / genotypesTMP.length;
        }

        varianceX = JSci.maths.ArrayMath.variance(snpmeancorrectedgenotypes);
        return new Pair<double[], Double>(snpmeancorrectedgenotypes, varianceX);
    }

    private double[] permuteGenotypes(double[] genotypes, Random rand) {
        ArrayList<Double> availableGenotypes = new ArrayList<Double>();
        for (int i = 0; i < genotypes.length; i++) {
            availableGenotypes.add(genotypes[i]);
        }
        Collections.shuffle(availableGenotypes, rand);
        return Primitives.toPrimitiveArr(availableGenotypes.toArray(new Double[0]));
    }
}
