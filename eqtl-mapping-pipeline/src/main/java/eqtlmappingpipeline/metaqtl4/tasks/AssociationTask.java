/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl4.tasks;

import eqtlmappingpipeline.metaqtl4.WorkPackage;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.Callable;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.TriTyperExpressionData;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDataset;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.Descriptives;

/**
 *
 * @author harmjan
 */
public class AssociationTask implements Callable<WorkPackage> {

    WorkPackage wp = null;
    short[][][] expressionToGenotypeCouplings; // [m_gg.length][m_settings.nrPermutationsFDR + 1][ind]
    private final int lengthOfCorrelationToZScoreTable;
    private final TriTyperGeneticalGenomicsDataset[] m_gg;
    private boolean cisOnly;
    private boolean transOnly;
    private boolean cisTrans;
    private final Integer[][] m_probeTranslation;
    private final int m_numProbes;

    public AssociationTask(WorkPackage wp,
            short[][][] expressionToGenotypeCouplings,
            int lengthOfCorrelationToZScoreTable,
            TriTyperGeneticalGenomicsDataset[] m_gg,
            boolean cis,
            boolean trans,
            Integer[][] probeTranslation,
            int numberTotalProbes) {
        this.wp = wp;
        this.expressionToGenotypeCouplings = expressionToGenotypeCouplings;
        this.lengthOfCorrelationToZScoreTable = lengthOfCorrelationToZScoreTable;
        this.m_gg = m_gg;
        cisOnly = false;
        cisTrans = false;
        transOnly = false;

        if (cis && !trans) {
            cisOnly = true;
        } else if (!cis && trans) {
            transOnly = true;
        } else if (cis && trans) {
            cisTrans = true;
        }

        this.m_probeTranslation = probeTranslation;
        this.m_numProbes = numberTotalProbes;
    }

    public WorkPackage call() throws Exception {
        // normalize genotype data
        // normalize once, test many
        SNP[] snps = wp.getSNPs();
        Triple<double[][], double[], boolean[][]> normalizedGenotypes = normalizeGenotypes(snps);
        double[][] snpgenotypes = normalizedGenotypes.getLeft();
        double[] snpvariances = normalizedGenotypes.getMiddle();
        boolean[][] includeExpressionSamples = normalizedGenotypes.getRight();
        int m_numDatasets = snps.length;
        int nrPermutations = expressionToGenotypeCouplings[0].length;

        TriTyperExpressionData[] m_expressiondata = new TriTyperExpressionData[m_numDatasets];
        for (int i = 0; i < m_expressiondata.length; i++) {
            m_expressiondata[i] = m_gg[i].getExpressionData();
        }

        int[] probes = wp.getProbes();
        HashSet<Integer> probestoExclude = null;
        if (transOnly) {
            if (probes != null) {
                probestoExclude = new HashSet<Integer>();
                for (int p = 0; p < probes.length; p++) {
                    probestoExclude.add(probes[p]);
                }
            }
        }

        Integer[][] pvalueIndexes = new Integer[nrPermutations+1][m_numProbes];
        for (int perm = 0; perm < nrPermutations; perm++) {
            double[][] zscores = null;
            int[] numSamples = new int[m_numDatasets];
            if (cisOnly) {
                zscores = new double[m_numDatasets][probes.length];
            } else if (transOnly) {
                zscores = new double[m_numDatasets][m_numProbes];;
            } else {
                zscores = new double[m_numDatasets][m_numProbes];
            }

            for (int d = 0; d < m_numDatasets; d++) {
                double[] snpmeancorrectedgenotypes = null;
                if (perm > 0) {
                    // permute the genotype data..
                    short[] realCouplings = expressionToGenotypeCouplings[d][0];
                    short[] newCouplings = expressionToGenotypeCouplings[d][perm]; // [m_gg.length][m_settings.nrPermutationsFDR + 1][ind]
                    double[] genotypes = snpgenotypes[d];
                    snpmeancorrectedgenotypes = new double[genotypes.length];

                    HashMap<Short, Short> genotypeToCtrId = new HashMap<Short, Short>();
                    short ctr = 0;
                    for (short i = 0; i < realCouplings.length; i++) {
                        if (includeExpressionSamples[d][i]) {
                            genotypeToCtrId.put(realCouplings[i], ctr);
                            ctr++;
                        }
                    }


                    // determine permuted expression to genotype coupling

                    // determine whether the expression sample was included for calculation

                    // if so, 

                    // determine position in original mappings

                    ctr = 0;
                    for (int i = 0; i < newCouplings.length; i++) {
                        short permId = newCouplings[i];
                        Short ctrId = genotypeToCtrId.get(permId);
                        if (ctrId != null) {
                            snpmeancorrectedgenotypes[ctr] = genotypes[ctrId];
                            ctr++;
                        }
                    }
                } else {
                    snpmeancorrectedgenotypes = snpgenotypes[d];
                }
//                snpmeancorrectedgenotypes = snpgenotypes[d];

                double snpvariance = snpvariances[d];
                boolean[] includeExpressionSample = includeExpressionSamples[d];
                double[] dsZScore = zscores[d];
                
                SNP dSNP = snps[d];
                if (dSNP != null) {
                    double[][] rawData = m_expressiondata[d].getMatrix();
                    double[] varY = m_expressiondata[d].getProbeVariance();
                    double[] meanY = m_expressiondata[d].getProbeMean();
                    int samplecount = m_expressiondata[d].getIndividuals().length;

                    numSamples[d] = snpmeancorrectedgenotypes.length;
                    if (cisOnly) {
                        for (int p = 0; p < probes.length; p++) {
                            int pid = probes[p];
                            Integer probeId = m_probeTranslation[d][pid];
                            if (probeId != null) {
                                double zScore = test(probeId, snpmeancorrectedgenotypes, snpvariance, varY[probeId], meanY[probeId], includeExpressionSample, samplecount, rawData);
                                dsZScore[p] = zScore;
                            }
                        }
                    } else if (transOnly) {
                        for (int pid = 0; pid < m_numProbes; pid++) {
                            if (probestoExclude == null || !probestoExclude.contains(pid)) {
                                Integer probeId = m_probeTranslation[d][pid];
                                if (probeId != null) {
                                    double zScore = test(probeId, snpmeancorrectedgenotypes, snpvariance, varY[probeId], meanY[probeId], includeExpressionSample, samplecount, rawData);
                                    dsZScore[pid] = zScore;
                                }
                            }
                        }
                    } else {
                        for (int pid = 0; pid < m_numProbes; pid++) {
                            Integer probeId = m_probeTranslation[d][pid];
                            if (probeId != null) {
                                double zScore = test(probeId, snpmeancorrectedgenotypes, snpvariance, varY[probeId], meanY[probeId], includeExpressionSample, samplecount, rawData);
                                dsZScore[pid] = zScore;
                            }
                        }
                    }
                }
            }

            if (perm == 0) {
                metaAnalyzeAndMergeData(pvalueIndexes[perm], numSamples, zscores);
                
            } else {
                metaAnalyzeAndMergeData(pvalueIndexes[perm], numSamples, zscores);
                
            }
        }
        clearSNPs();
        wp.setResults(pvalueIndexes);
        return wp;
    }

    private double test(Integer probeId, double[] x, double varianceX, double varianceY, double meanY, boolean[] includeExpressionSample, int sampleCount, double[][] rawData) {
        double[] y = null;

//        long t1 = System.nanoTime();
//        long t2 = 0;
        if (x.length != sampleCount) {
            y = new double[x.length];
            int itr = 0;
            double sum = 0;
            double[] tmpY = rawData[probeId];
            // recalculate mean and variance
            for (int s = 0; s < sampleCount; s++) {
                if (includeExpressionSample[s]) {
                    y[itr] = tmpY[s];
                    sum += y[itr];
                    itr++;
                }
            }

            meanY = sum / itr;

            double varsum = 0;
            for (int i = 0; i < y.length; i++) {
                y[i] -= meanY;
                varsum += y[i] * y[i];
            }
            varianceY = varsum / (y.length - 1);
//            t2 = System.nanoTime();
        } else {
//            y = rawData[probeId];
            y = new double[x.length];
            System.arraycopy(rawData[probeId], 0, y, 0, x.length);
        }

        //Calculate correlation coefficient:
        if (varianceY == 0) {
            return Double.NaN;
        } else {
            double correlation = Correlation.correlate(x, y, varianceX, varianceY);
            if (correlation >= -1 && correlation <= 1) {
                double zScore = Correlation.convertCorrelationToZScore(x.length, correlation);
                return zScore;
            } else {
                System.err.println("Error! correlation invalid: " + correlation);
                System.exit(-1);
            }
        }
        return Double.NaN;
    }

    private void clearSNPs() {
        SNP[] snps = wp.getSNPs();
        for (SNP s : snps) {
            if (s != null) {
                s.clearGenotypes();
            }
        }
    }

    private Triple<double[][], double[], boolean[][]> normalizeGenotypes(SNP[] snps) {

        int m_numDatasets = snps.length;
        double[] snpvariances = new double[m_numDatasets];
        double[][] snpmeancorrectedgenotypes = new double[m_numDatasets][0];
        boolean[][] includeExpressionSample = new boolean[m_numDatasets][0];

        for (int d = 0; d < m_numDatasets; d++) {
            SNP dSNP = snps[d];

            if (dSNP != null) {

                short[] expressionToGenotypeId = expressionToGenotypeCouplings[d][0]; // [m_gg.length][m_settings.nrPermutationsFDR + 1][ind]
                double[] x = dSNP.selectGenotypes(expressionToGenotypeId, false, true);

                int xLen = x.length;
                double meanX = JSci.maths.ArrayMath.mean(x);

                snpmeancorrectedgenotypes[d] = new double[xLen];

                for (int i = 0; i < xLen; i++) {
                    snpmeancorrectedgenotypes[d][i] = x[i] - meanX;
                }

                double varianceX = JSci.maths.ArrayMath.variance(x);
                if (varianceX != 0) {
                    snpvariances[d] = varianceX;

                    short inds[] = expressionToGenotypeId;
                    int sampleCount = expressionToGenotypeId.length;
                    includeExpressionSample[d] = new boolean[sampleCount];
                    byte[] genotypes = dSNP.getGenotypes();
                    for (int s = 0; s < sampleCount; s++) {
                        int ind = inds[s];
                        double valX = genotypes[ind]; // loadedSNPGenotype[ind];
                        if (valX != -1) {
                            includeExpressionSample[d][s] = true;
                        } else {
                            includeExpressionSample[d][s] = false;
                        }
                    }
                } else {
                    dSNP.clearGenotypes();
                    dSNP = null;
                    wp.getFlipSNPAlleles()[d] = null;
                    snps[d] = null;
                }
            }
        }
        return new Triple<double[][], double[], boolean[][]>(snpmeancorrectedgenotypes, snpvariances, includeExpressionSample);
    }

    private void metaAnalyzeAndMergeData(Integer[] pvaluedist, int[] numSamplesPerDataset, double[][] zscores) {
        // meta analyze z-scores
        // per probe, convert to p-value
        int numProbes = zscores[0].length;
        int m_numDatasets = wp.getSNPs().length;

        for (int p = 0; p < numProbes; p++) {
            int nrDatasetsPassingQC = 0;
            int nrTotalSamples = 0;
            double zSum = 0;

            for (int d = 0; d < m_numDatasets; d++) {
                double zscore = zscores[d][p];
                if (!Double.isNaN(zscore)) {
                    int numSamples = numSamplesPerDataset[d];
                    boolean flipalleles = wp.getFlipSNPAlleles()[d];
                    if (flipalleles) {
                        zscore = -zscore;
                    }
                    nrDatasetsPassingQC++;

                    double weight = Descriptives.getSqrt(numSamples);

                    zSum += (zscore * weight);
                    nrTotalSamples += numSamples;
                }
            }
            if (nrDatasetsPassingQC > 0 && nrTotalSamples > 0) {
                double sqrtSample = Descriptives.getSqrt(nrTotalSamples);
                double zScore = zSum / sqrtSample;

                int index = Descriptives.getZScorePvalueIndex(zScore);
                int probeId = p;
                if(cisOnly){
                    probeId = wp.getProbes()[p];
                }
                pvaluedist[probeId] = index;
            }
        }
    }
}
