/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl4;

import gnu.trove.map.hash.TObjectIntHashMap;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.Callable;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

/**
 *
 * @author harm-jan
 */
public class MetaQTL4CorrelationTask implements Callable<Pair<int[], int[]>> {

    private final int distributionSize;
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
    private final int threadIndex;
    private final int[] realDistribution;
    private final int[] permutedDistribution;

    MetaQTL4CorrelationTask(int nrThreads,
            int distributionSize,
            long[] randomizationSeeds,
            ArrayList<MetaQTL4MetaTrait> availableTraits,
            TObjectIntHashMap<MetaQTL4MetaTrait> availableTraitsHash,
            MetaQTL4Dataset[] datasets,
            GeneticVariant[][] geneticVariantIndex,
            MetaQTL4Settings m_settings,
            MetaQTL4TraitAnnotation traitAnnotation,
            Integer[][] traitIndex,
            Set<MetaQTL4MetaTrait> traitsToInclude,
            Set<String> variantsToInclude,
            int threadIndex) {
        this.distributionSize = distributionSize;
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
        this.threadIndex = threadIndex;
        this.realDistribution = new int[distributionSize];
        this.permutedDistribution = new int[distributionSize];
    }

    @Override
    public Pair<int[], int[]> call() throws Exception {
        try {
            // normal procedure
            int nrPermutations = randomizationSeeds.length;
            int nrDatasets = geneticVariantIndex.length;
            int nrVariants = geneticVariantIndex[nrDatasets - 1].length;

            // iterate the variants
            ArrayList<Integer> traitsToTest = null;
            boolean[] flipEffects = new boolean[nrDatasets];
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
            for (int variantId = 0; variantId < nrVariants; variantId += threadIndex) {
                pb.iterate();
                // load the genotypes
                double[][] genotypes = new double[nrDatasets][0];
                double[] genotypeVariances = new double[nrDatasets];
                int[] sampleSizes = new int[nrDatasets];

                int startPosition = -1;
                String sequence = null;
                numberOfDatasetsPassingQC = 0;

                for (int datasetId = 0; datasetId < nrDatasets; datasetId++) {
                    GeneticVariant variant = geneticVariantIndex[datasetId][variantId];

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
                            genotypes[datasetId] = genotypedata.getLeft();

                            genotypeVariances[datasetId] = genotypedata.getRight();
                            sampleSizes[datasetId] = genotypes[datasetId].length;
                            if (sequence == null) {
                                startPosition = variant.getStartPos();
                                sequence = variant.getSequenceName();
                            }
                        } else {
                            // TODO: output to the log?
                            genotypes[datasetId] = null;
                        }
                    } else {
                        genotypes[datasetId] = null;
                    }
                }

                // determine which traits to test
                // If this is a cis OR a trans analysis but NOT BOTH: use XOR
                // else test ALL traits
                if (m_settings.isCisAnalysis() ^ m_settings.isTransAnalysis()) {
                    // index once, use many
                    if (traitsToTest == null) {
                        traitsToTest = new ArrayList<Integer>();
                        // possible faster with hashmap lookup
                        for (int i = 0; i < availableTraits.size(); i++) {
                            MetaQTL4MetaTrait trait = availableTraits.get(i); // a treemap may be more efficient here
                            boolean sameChr = trait.getChr().equals(sequence);
                            if (sameChr) {
                                int distance = Math.abs(trait.getChrMidpoint() - startPosition); // precalculate this, if possible
                                if (m_settings.isCisAnalysis() && distance < m_settings.getCiseQTLAnalysMaxSNPProbeMidPointDistance()) {
                                    traitsToTest.add(i);
                                }
                                if (m_settings.isTransAnalysis() && distance > m_settings.getCiseQTLAnalysMaxSNPProbeMidPointDistance()) {
                                    traitsToTest.add(i);
                                }
                            }
                            if (!sameChr && m_settings.isTransAnalysis()) {
                                traitsToTest.add(i);
                            }
                        }
                    }
                }

                if (traitsToTest != null && traitsToTest.isEmpty()) {
                    // do something, log perhaps
                } else {
                    // run test
                    int[] dist;
                    for (int permutation = -1; permutation < nrPermutations; permutation++) {

                        if (permutation >= 0) {
                            dist = permutedDistribution;
                            long seed = randomizationSeeds[permutation];
                            Random rand = new Random(seed);

                            // reorder genotypes on the basis of the permuted sample links
                            for (int d = 0; d < datasets.length; d++) {
                                if (genotypes[d] != null) {
                                    genotypes[d] = permuteGenotypes(genotypes[d], rand);
                                }
                            }
                        } else {
                            dist = realDistribution;
                        }

                        if (traitsToTest == null) {
                            // cistrans analysis
                            for (int trait = 0; trait < availableTraits.size(); trait++) {
                                int bin = test(sampleSizes, genotypes, genotypeVariances, trait, includeTraitSample);
                                dist[bin]++;
                            }
                        } else {
                            for (Integer trait : traitsToTest) {
//                                MetaQTL4MetaTrait traitObj = availableTraits.get(trait);
//                                System.out.println(Strings.concat(traitObj.getPlatformIds(), Strings.comma));
                                int bin = test(sampleSizes, genotypes, genotypeVariances, trait, includeTraitSample);
                                dist[bin]++;
                            }
                        }
                    }
                }
            }
            pb.close();

        } catch (Exception e) {
            e.printStackTrace();
        }
        return new Pair<int[], int[]>(realDistribution, permutedDistribution);
    }

    private int test(int[] sampleSizes, double[][] genotypes, double[] genotypeVariances, Integer traitId, boolean[][] includeTraitSample) {
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
                    if (sampleSizes[datasetId] != y.length) {
                        double[] newY = new double[x.length];
                        int itr = 0;
                        double sum = 0;

                        // recalculate mean and variance
                        for (int s = 0; s < y.length; s++) {
                            if (includeDatasetTraitSample[s]) {
                                newY[itr] = y[s];
                                sum += newY[itr];
                                itr++;
                            }
                        }

                        y = newY;
                        meanY = sum / itr;
                        double varsum = 0;
                        for (int i = 0; i < y.length; i++) {
                            y[i] -= meanY;
                            varsum += y[i] * y[i];
                        }
                        varianceY = varsum / (y.length - 1);
                    } else {
                        varianceY = datasets[datasetId].getTraitVariance(datasetTraitId);
                    }

                    if (varianceY == 0) {
                        // trait has no variance
                        zscores[datasetId] = Double.NaN;
                    } else {
                        //Calculate correlation coefficient:
                        double correlation = Correlation.correlate(x, y, varianceX, varianceY);
//                        double correlation2 = JSci.maths.ArrayMath.correlation(x, y);
//
//                        if (correlation != correlation2) {
//                            System.err.println("ERROR: "+(correlation-correlation2));
//                            double[] newY = new double[x.length];
//                            int itr = 0;
//                            double sum = 0;
//
//                            // recalculate mean and variance
//                            for (int s = 0; s < y.length; s++) {
//                                if (includeDatasetTraitSample[s]) {
//                                    newY[itr] = y[s];
//                                    sum += newY[itr];
//                                    itr++;
//                                }
//                            }
//
//                            meanY = sum / itr;
//                            double varsum = 0;
//                            for (int i = 0; i < newY.length; i++) {
//                                newY[i] -= meanY;
//                                varsum += newY[i] * newY[i];
//                            }
//                            double varianceYRecalculated = varsum / (newY.length - 1);
//                            System.err.println(correlation + "\t" + correlation2 + "\t"
//                                    + x.length + "\t"
//                                    + varianceYRecalculated + "\t"
//                                    + varianceX + "\t"
//                                    + varianceY + "\t"
//                                    + datasets[datasetId].getTraitVariance(datasetTraitId) + "\t"
//                                    + Descriptives.variance(x) + "\t"
//                                    + Descriptives.variance(y)
//                            );
////                            for (int i = 0; i < y.length; i++) {
////                                System.out.println(i + "\t" + y[i] + "\t" + x[i]);
////                            }
////                            System.exit(-1);
//                        }

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
        int bin = (int) Math.round(p * distributionSize);
        if (bin == distributionSize) {
            bin--;
        }

        return bin;
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

        if (genotypes == null) {
            System.err.println("ERROR: genotypes null");
            return null;
        }

        for (int i = 0; i < genotypes.length; i++) {
            availableGenotypes.add(genotypes[i]);
        }
        Collections.shuffle(availableGenotypes, rand);
        return Primitives.toPrimitiveArr(availableGenotypes.toArray(new Double[0]));
    }

}
