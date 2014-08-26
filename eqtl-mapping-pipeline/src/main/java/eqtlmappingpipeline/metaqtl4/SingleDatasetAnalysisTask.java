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
import umcg.genetica.containers.Pair;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.util.Primitives;

/**
 *
 * @author Harm-Jan
 */
public class SingleDatasetAnalysisTask implements Callable<Boolean> {

    private final long[] randomizationSeeds;
    private final ArrayList<MetaQTL4MetaTrait> availableTraits;
    private final MetaQTL4MetaTraitTreeSet availableTraitsHash;
    private final MetaQTL4Dataset dataset;
    private final MetaQTL4Settings m_settings;

    private final CompletionService resultPool;
    private final GeneticVariant variant;

    double mafthreshold;
    double hwepthreshold;
    double callratethreshold;
    Integer[] traitIndex;

    public SingleDatasetAnalysisTask(int nrThreads, long[] randomizationSeeds,
            ArrayList<MetaQTL4MetaTrait> availableTraits, MetaQTL4MetaTraitTreeSet availableTraitsHash, Integer[] traitIndex,
            MetaQTL4Dataset dataset, MetaQTL4Settings m_settings, Set<MetaQTL4MetaTrait> traitsToInclude,
            GeneticVariant variant, CompletionService resultPool) {

        this.randomizationSeeds = randomizationSeeds;
        this.availableTraits = availableTraits;
        this.availableTraitsHash = availableTraitsHash;
        this.dataset = dataset;
        this.m_settings = m_settings;

        this.traitIndex = traitIndex;
        this.resultPool = resultPool;
        this.variant = variant;

        mafthreshold = m_settings.getSnpQCMAFThreshold();
        hwepthreshold = m_settings.getSnpQCHWEThreshold();
        callratethreshold = m_settings.getSnpQCHWEThreshold();

    }

    @Override
    public Boolean call() throws Exception {
        try {
//            System.out.println("Running " + variant.getPrimaryVariantId());
            // normal procedure
            int nrPermutations = randomizationSeeds.length;

            // iterate the variants
            int maxNrSamples = dataset.getGenotypeToTraitCouplingInt().length;
            boolean[] includeTraitSample = new boolean[maxNrSamples];
            float[] genotypesTMP = dataset.getSampleDosages(variant);
            if (!variant.isBiallelic() || variant.getHwePvalue() < hwepthreshold || variant.getMinorAlleleFrequency() < mafthreshold || variant.getCallRate() < callratethreshold) {
                return true;
            }
            Pair<double[], Double> genotypedata = correctGenotypesForMissingValuesAndNormalize(dataset.getGenotypeToTraitCouplingInt(), variant, genotypesTMP, includeTraitSample);
            int snpStartPosition = variant.getStartPos();
            String sequence = variant.getSequenceName();
            // determine list of probes for each variant
            HashSet<Integer> traitsToTest = null; //

            double[] datasetGenotypeData = genotypedata.getLeft();
            int sampleSize = datasetGenotypeData.length;
            double genotypeVariance = genotypedata.getRight();

            if (m_settings.isCisAnalysis() ^ m_settings.isTransAnalysis()) {
                // index once, use many
                traitsToTest = new HashSet<Integer>();

//                for (int i = 0; i < availableTraits.size(); i++) {
//                    MetaQTL4MetaTrait trait = availableTraits.get(i); // a treemap may be more efficient here
//                    boolean sameChr = trait.getChr().equals(sequence);
//                    if (sameChr) {
//                        int distance = Math.abs(trait.getChrMidpoint() - snpStartPosition); // precalculate this, if possible
//                        if (m_settings.isCisAnalysis() && distance < m_settings.getCiseQTLAnalysMaxSNPProbeMidPointDistance()) {
//                            traitsToTest.add(i);
//                        }
//                        if (m_settings.isTransAnalysis() && distance < m_settings.getCiseQTLAnalysMaxSNPProbeMidPointDistance()) {
//                            traitsToTest.add(i);
//                        }
//                    }
//                }
                Set<MetaQTL4MetaTrait> traits = availableTraitsHash.getTraitInWindow(sequence, snpStartPosition, m_settings.getCiseQTLAnalysMaxSNPProbeMidPointDistance()); // a treemap may be more efficient here
                for (MetaQTL4MetaTrait trait : traits) {
                    traitsToTest.add(trait.getCurrentMetaId());
                }
            }

            // run the analysis
            for (int permutation = -1; permutation < nrPermutations; permutation++) {
                // iterate the variants per dataset..

                // permute the genotypes, if required.
                if (permutation >= 0) {
                    long seed = randomizationSeeds[permutation];
                    Random rand = new Random(seed);
                    // reorder genotypes on the basis of the permuted sample links

                    if (datasetGenotypeData != null) {
                        datasetGenotypeData = permuteGenotypes(datasetGenotypeData, rand);
                    }

                }

                if (traitsToTest == null) {
                    // cistrans analysis
                    for (int trait = 0; trait < availableTraits.size(); trait++) {
                        //int bin = test(sampleSizes, genotypes, genotypeVariances, trait, includeTraitSample);
                        test(sampleSize, datasetGenotypeData, genotypeVariance, trait, includeTraitSample);
                    }
                } else {

                    if (m_settings.isCisAnalysis()) {

                        for (Integer trait : traitsToTest) {
                            test(sampleSize, datasetGenotypeData, genotypeVariance, trait, includeTraitSample);
                        }
                    }
                    if (m_settings.isTransAnalysis()) {
                        for (int i = 0; i < availableTraits.size(); i++) {
                            if (!traitsToTest.contains(i)) {
                                test(sampleSize, datasetGenotypeData, genotypeVariance, i, includeTraitSample);
                            }
                        }

                    }

                }

            }

        } catch (Exception e) {
            e.printStackTrace();
        }
        return true;
    }

    private void test(int sampleSize, double[] genotypes, double genotypeVariance, Integer traitId, boolean[] includeTraitSample) {

        double zscore = 0;

        if (genotypes == null) {
            zscore = Double.NaN;
        } else {
            double meanY;
            double varianceY;

            // re-normalize the trait data when the genotypes had missing values
            Integer datasetId = traitIndex[traitId];
            double[] y = dataset.getTraitData(datasetId);
//            if (sampleSize != y.length) {
            double[] newY = new double[genotypes.length];
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
            double varianceYOrig = dataset.getTraitVariance(traitId);

            if (varianceY == 0) {
                // trait has no variance
                zscore = Double.NaN;
            } else {
                //Calculate correlation coefficient:
                double correlation = Correlation.correlateMeanCenteredData(genotypes, y, genotypeVariance, varianceY);
                //double correlation2 = JSci.maths.ArrayMath.correlation(genotypes, y);
                //System.out.println(correlation + "\t" + correlation2 + "\t" + genotypes.length + "\t" + genotypeVariance + "\t" + varianceY + "\t" + varianceY + "\t" + Descriptives.variance(genotypes) + "\t" + Descriptives.variance(y));
                if (correlation >= -1 && correlation <= 1) {
                    zscore = Correlation.convertCorrelationToZScore(genotypes.length, correlation);
                    
                    double p = Descriptives.convertZscoreToPvalue(zscore);
                } else {
                    System.err.println("Error! correlation invalid: " + correlation);
                    for (int i = 0; i < genotypes.length; i++) {
                        System.out.println(genotypes[i] + "\t" + y[i] + "\t" + newY[i]);
                    }
                    System.exit(-1);
                }

            }

        }

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
