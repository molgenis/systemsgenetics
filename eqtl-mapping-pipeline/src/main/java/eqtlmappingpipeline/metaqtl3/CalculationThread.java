/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl3;

import eqtlmappingpipeline.metaqtl3.containers.WorkPackage;
import eqtlmappingpipeline.metaqtl3.containers.Result;
import umcg.genetica.math.stats.Descriptives;
import eqtlmappingpipeline.metaqtl3.graphics.EQTLPlotter;
import java.nio.ByteBuffer;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.zip.Deflater;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.TriTyperExpressionData;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.Regression;
import umcg.genetica.util.RunTimer;

/**
 *
 * @author harmjan
 */
class CalculationThread extends Thread {

    final TriTyperExpressionData[] m_expressiondata;
    final Integer[][] m_probeTranslation;
    int m_name;
    private int m_numProbes;
    private int m_numDatasets;
    private final int[][] m_expressionToGenotypeIds;
    private final double[][] probeVariance;
    private final double[][] probeMean;
    private final String[][] probeName;
    private final LinkedBlockingQueue<WorkPackage> m_workpackage_queue;
    private final LinkedBlockingQueue<WorkPackage> m_result_queue;
    int testsPerformed = 0;
    public boolean done = false;
    private int failedQC;
    private boolean cisOnly;
    private boolean cisTrans;
    private boolean transOnly;
    private boolean useAbsolutePValues;
    private final EQTLPlotter m_eQTLPlotter;
    private final double m_pvaluePlotThreshold;
    private boolean determinebeta = false;
    private boolean determinefoldchange = false;
    private WorkPackage currentWP;
    private boolean m_binaryoutput = false;

    CalculationThread(int i, LinkedBlockingQueue<WorkPackage> packageQueue, LinkedBlockingQueue<WorkPackage> resultQueue, TriTyperExpressionData[] expressiondata, Integer[][] probeTranslationTable,
            int[][] expressionToGenotypeIds, MetaQTL3Settings settings, EQTLPlotter plotter, boolean binaryoutput) {
        m_binaryoutput = binaryoutput;
        m_name = i;
        m_workpackage_queue = packageQueue;
        m_result_queue = resultQueue;
        m_probeTranslation = probeTranslationTable;
        m_expressiondata = expressiondata;
        boolean m_cis = settings.cisAnalysis;
        boolean m_trans = settings.transAnalysis;
        m_name = i;
        m_numProbes = m_probeTranslation[m_probeTranslation.length - 1].length;
        m_numDatasets = m_probeTranslation.length;
        m_expressionToGenotypeIds = expressionToGenotypeIds;
        probeVariance = new double[m_numDatasets][0];
        probeMean = new double[m_numDatasets][0];
        probeName = new String[m_numDatasets][0];
        for (int d = 0; d < m_numDatasets; d++) {
            probeVariance[d] = expressiondata[d].getProbeVariance();
            probeMean[d] = expressiondata[d].getProbeMean();
            probeName[d] = expressiondata[d].getProbes();
        }

        cisOnly = false;
        cisTrans = false;
        transOnly = false;

        determinebeta = settings.provideBetasAndStandardErrors;
        determinefoldchange = settings.provideFoldChangeData;

        if (m_cis && !m_trans) {
            cisOnly = true;
        } else if (!m_cis && m_trans) {
            transOnly = true;
        } else if (m_cis && m_trans) {
            cisTrans = true;
        }

        m_eQTLPlotter = plotter;
        m_pvaluePlotThreshold = settings.plotOutputPValueCutOff;
    }

    @Override
    public void run() {
//        System.out.println("Starting thread: " + this.getName() +". Waiting for input");

        // while the process has not been killed or when the queue is not empty
        boolean poison = false;

        RunTimer t = new RunTimer();
        t.start();
        int taken = 0;
        while (!poison) {
            try {
                WorkPackage pack = m_workpackage_queue.take();
                if (!pack.getPoison()) {
                    analyze(pack);
                    taken++;

//                    if(taken % printperiterations == 0){
//                        System.out.println("Thread "+this.getName()+" calculated "+taken+" workpackages.");
//                    }
                } else {
                    poison = pack.getPoison();
//                    System.out.println("Thread " + m_name + " got killed by a poisonous workpackage, but was bravely able to perform\t" + testsPerformed + "\ttests");
                }

            } catch (InterruptedException ex) {
                ex.printStackTrace();
            }
        }
//        System.out.println(this.getName() + " finished.\t" + testsPerformed + "\ttests performed");
    }

    public void kill() {
        done = true;
    }

    private void analyze(WorkPackage wp) {
        testsPerformed = 0;
        currentWP = wp;
        wp.setNumTested(0);
//        RunTimer t1 = new RunTimer();
        // load SNP genotypes
        SNP[] snps = wp.getSnps();
        int[] probes = wp.getProbes();
        Result dsResults = null;

        double[] snpvariances = new double[m_numDatasets];
        double[][] snpmeancorrectedgenotypes = new double[m_numDatasets][0];
        double[][] originalgenotypes = new double[m_numDatasets][0];
        boolean[][] includeExpressionSample = new boolean[m_numDatasets][0];

        for (int d = 0; d < m_numDatasets; d++) {
            SNP dSNP = snps[d];

            if (dSNP != null) {

                double[] x = dSNP.selectGenotypes(m_expressionToGenotypeIds[d], false, true);
                originalgenotypes[d] = dSNP.selectGenotypes(m_expressionToGenotypeIds[d], false, false);

                int xLen = x.length;
                double meanX = JSci.maths.ArrayMath.mean(x);

                snpmeancorrectedgenotypes[d] = new double[xLen];

                for (int i = 0; i < xLen; i++) {
                    snpmeancorrectedgenotypes[d][i] = x[i] - meanX;
                }

                double varianceX = JSci.maths.ArrayMath.variance(x);
                if (varianceX != 0) {
                    snpvariances[d] = varianceX;

                    int inds[] = m_expressionToGenotypeIds[d];
                    int sampleCount = m_expressionToGenotypeIds[d].length;
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

        if (cisOnly) {
            dsResults = new Result(m_numDatasets, wp.getProbes().length, wp.getId());
            for (int d = 0; d < m_numDatasets; d++) {
                SNP dSNP = snps[d];

                if (dSNP != null) {
                    dsResults.numSamples[d] = snpmeancorrectedgenotypes[d].length;
                    double[][] rawData = m_expressiondata[d].getMatrix();
                    double[] varY = m_expressiondata[d].getProbeVariance();
                    double[] meanY = m_expressiondata[d].getProbeMean();
                    int samplecount = m_expressiondata[d].getIndividuals().length;
                    for (int p = 0; p < probes.length; p++) {
                        int pid = probes[p];
                        Integer probeId = m_probeTranslation[d][pid];
                        if (probeId != null) {
                            test(d, p, probeId, snpmeancorrectedgenotypes[d], originalgenotypes[d], snpvariances[d], varY[probeId], meanY[probeId], includeExpressionSample[d], samplecount, rawData, dsResults);
                        } else {
                            dsResults.correlations[d][p] = Double.NaN;
                            dsResults.zscores[d][p] = Double.NaN;
                        }
                    }

                } else {
                    for (int p = 0; p < probes.length; p++) {
                        dsResults.correlations[d][p] = Double.NaN;
                        dsResults.zscores[d][p] = Double.NaN;
                    }
                }
            }
        } else if (transOnly) {

            HashSet<Integer> probestoExclude = null;
            if (probes != null) {
                probestoExclude = new HashSet<Integer>();
                for (int p = 0; p < probes.length; p++) {
                    probestoExclude.add(probes[p]);
                }
            }
            dsResults = new Result(m_numDatasets, m_numProbes, wp.getId());
            for (int d = 0; d < m_numDatasets; d++) {
                SNP dSNP = snps[d];
                dsResults.numSamples[d] = snpmeancorrectedgenotypes[d].length;
                double[][] rawData = m_expressiondata[d].getMatrix();
                double[] varY = m_expressiondata[d].getProbeVariance();
                double[] meanY = m_expressiondata[d].getProbeMean();
                int samplecount = m_expressiondata[d].getIndividuals().length;
                if (dSNP != null) {
                    dsResults.numSamples[d] = snpmeancorrectedgenotypes[d].length;
                    for (int pid = 0; pid < m_numProbes; pid++) {
                        if (probestoExclude == null || !probestoExclude.contains(pid)) {
                            Integer probeId = m_probeTranslation[d][pid];
                            if (probeId != null) {
                                test(d, pid, probeId, snpmeancorrectedgenotypes[d], originalgenotypes[d], snpvariances[d], varY[probeId], meanY[probeId], includeExpressionSample[d], samplecount, rawData, dsResults);
                            } else {
                                dsResults.correlations[d][pid] = Double.NaN;
                                dsResults.zscores[d][pid] = Double.NaN;
                            }
                        } else {
                            dsResults.correlations[d][pid] = Double.NaN;
                            dsResults.zscores[d][pid] = Double.NaN;
                        }
                    }
                } else {
                    for (int p = 0; p < m_numProbes; p++) {
                        dsResults.correlations[d][p] = Double.NaN;
                        dsResults.zscores[d][p] = Double.NaN;
                    }
                }

            }
        } else {
            dsResults = new Result(m_numDatasets, m_numProbes, wp.getId());
            for (int d = 0; d < m_numDatasets; d++) {
                SNP dSNP = snps[d];
                dsResults.numSamples[d] = snpmeancorrectedgenotypes[d].length;
                double[][] rawData = m_expressiondata[d].getMatrix();
                double[] varY = m_expressiondata[d].getProbeVariance();
                double[] meanY = m_expressiondata[d].getProbeMean();
                int samplecount = m_expressiondata[d].getIndividuals().length;
                if (dSNP != null) {
                    dsResults.numSamples[d] = snpmeancorrectedgenotypes[d].length;
//                    RunTimer t2 = new RunTimer();
                    for (int pid = 0; pid < m_numProbes; pid++) {
                        Integer probeId = m_probeTranslation[d][pid];
                        if (probeId != null) {
                            test(d, pid, probeId, snpmeancorrectedgenotypes[d], originalgenotypes[d], snpvariances[d], varY[probeId], meanY[probeId], includeExpressionSample[d], samplecount, rawData, dsResults);
                        } else {
                            dsResults.correlations[d][pid] = Double.NaN;
                            dsResults.zscores[d][pid] = Double.NaN;
                        }
                    }
//                    System.out.println("Test: "+t2.getTimeDesc());
                } else {
                    for (int p = 0; p < m_numProbes; p++) {
                        dsResults.correlations[d][p] = Double.NaN;
                        dsResults.zscores[d][p] = Double.NaN;
                    }
                }

            }
        }

        convertResultsToPValues(wp, dsResults);

        if (m_eQTLPlotter != null) {
            for (int p = 0; p < dsResults.pvalues.length; p++) {
                Double pval = dsResults.pvalues[p];
                if (pval != null && pval < m_pvaluePlotThreshold) {
                    ploteQTL(wp, p);
                }
            }
        }

        snps = wp.getSnps();
        if (snps != null) {
            for (int d = 0; d < snps.length; d++) {
                if (snps[d] != null) {
                    snps[d].clearGenotypes();
                }
            }
        }

        // if result output is binary, convert to bytes and deflate the set of bytes.

        if (m_binaryoutput) {
            deflateResults(wp);
        }
        // now push the results in the queue..
        try {
            wp.setNumTested(testsPerformed);
            m_result_queue.put(wp);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
//        System.out.println("Analyze: "+t1.getTimeDesc());
    }

    private void test(int d, int p, Integer probeId, double[] x, double[] originalGenotypes, double varianceX, double varianceY, double meanY, boolean[] includeExpressionSample, int sampleCount, double[][] rawData, Result r) {
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
            r.zscores[d][p] = Double.NaN;
            r.correlations[d][p] = Double.NaN;
        } else {
            //            long t3 = System.nanoTime();
            double correlation = Correlation.correlate(x, y, varianceX, varianceY);
            if (correlation >= -1 && correlation <= 1) {
                double zScore = Correlation.convertCorrelationToZScore(x.length, correlation);

//                if (determinebeta || determinefoldchange) {
                double[] xcopy = new double[x.length];

//                System.out.println("");
//                System.out.println("Pre");
//                for (int i = 0; i < y.length; i++) {
//                    System.out.println(i + "\t" + x[i] + "\t" + y[i]);
//                }

                double meanx = JSci.maths.ArrayMath.mean(x);
                double meany = JSci.maths.ArrayMath.mean(y);
                double stdevy = JSci.maths.ArrayMath.standardDeviation(y);
                double stdevx = JSci.maths.ArrayMath.standardDeviation(x);
                //
                for (int i = 0; i < y.length; i++) {
                    y[i] -= meany;
                    y[i] /= stdevy;
                    xcopy[i] = x[i] - meanx;
                    xcopy[i] /= stdevx;
                }
                meany = JSci.maths.ArrayMath.mean(y);
                meanx = JSci.maths.ArrayMath.mean(xcopy);
//                if (determinebeta) {
                calculateRegressionCoefficients(xcopy, meanx, y, meany, r, d, p);
//                }
                if (determinefoldchange) {
                    determineFoldchange(originalGenotypes, y, r, d, p);
                }
//                }

                //            long t4 = System.nanoTime();
                r.zscores[d][p] = zScore;
                r.correlations[d][p] = correlation;

//                System.out.println("");
//                System.out.println("Post");
//                for (int i = 0; i < y.length; i++) {
//                    System.out.println(i + "\t" + xcopy[i] + "\t" + y[i]);
//                }

//                double[] reg = Regression.getLinearRegressionCoefficients(xcopy, y);

//                System.out.println("z " + zScore + "\tcorr " + correlation + "\tbeta " + r.beta[d][p] + "\tse " + r.se[d][p] + "\t" + reg[0] + "\t" + reg[1] + "\t" + reg[2]);
//                System.out.println("");
//                System.exit(0);
            } else {
                System.err.println("Error! correlation invalid: " + correlation);
                System.exit(-1);
            }

//            long t5 = System.nanoTime();
//            System.out.println((t2-t1)+"\t"+(t3-t1)+"\t"+(t4-t1)+"\t"+(t5-t1));
        }
    }

    private void normalizeGeneExpressionDataForMissingGenotypes() {
    }

    private void calculateRegressionCoefficients(double[] x, double meanx, double[] y, double meany, Result r, int d, int p) {
        double beta = 0;
        double alpha = 0;
        double sxx = 0;
        double sxy = 0;
        double b = 0;

        for (int i = 0; i < y.length; i++) {
            sxx += ((x[i] - meanx) * (x[i] - meanx));
            sxy += ((y[i] - meany) * (x[i] - meanx));
        }

        beta = sxy / sxx;
        alpha = meany - beta * meanx;

        double ssxy = 0;
        for (int i = 0; i < y.length; i++) {
            double yexp = alpha + (beta * x[i]);
            ssxy += ((y[i] - yexp) * (y[i] - yexp));
        }

        double se = (Math.sqrt((ssxy) / (y.length - 2))) / Math.sqrt(sxx);
//	double se3 = se2 / Math.sqrt(sxx);

//            System.out.println("---");
//            System.out.println(beta+"\t"+alpha+"\t"+se2+"\t"+se3);
        r.beta[d][p] = beta;
        r.se[d][p] = se;
//            System.exit(0);

    }

    private void determineFoldchange(double[] genotypes, double[] expression, Result r, int d, int p) {
        int numAA = 0;
        int numBB = 0;

        double sumAA = 0;
        double sumBB = 0;

        for (int i = 0; i < genotypes.length; i++) {
//	    if (indWGA[i] != -1) {
            if (genotypes[i] == 0) {
                sumAA += (expression[i] * 2);
                numAA += 2;
            }
            if (genotypes[i] == 1) {
                sumAA += (expression[i]);
                sumBB += (expression[i]);
                numAA++;
                numBB++;
            }
            if (genotypes[i] == 2) {
                sumBB += (expression[i] * 2);
                numBB += 2;
            }
//	    }

        }





        sumAA /= (double) numAA;
        sumBB /= (double) numBB;

        double min = sumAA;
        if (sumBB < min) {
            min = sumBB;
        }

        // normalize to mean of 1
        if (min < 0) {
            sumAA += Math.abs(min) + 1;
            sumBB += Math.abs(min) + 1;
        }
        if (currentWP.getFlipSNPAlleles()[d]) {


            r.fc[d][p] = sumAA / sumBB;
        } else {
            r.fc[d][p] = sumBB / sumAA;
        }

    }

    private void convertResultsToPValues(WorkPackage wp, Result dsResults) {
        // per probe, convert to p-value
        int numProbes = dsResults.zscores[0].length;


        boolean hasResults = false;


        for (int p = 0; p < numProbes; p++) {
            int nrDatasetsPassingQC = 0;
            int nrTotalSamples = 0;
            double zSum = 0;
            double zSumAbsolute = 0;
            double betasum = 0;

            for (int d = 0; d < m_numDatasets; d++) {
                double zscore = dsResults.zscores[d][p];
                double correlation = dsResults.correlations[d][p];

                Integer numSamples = dsResults.numSamples[d];
                if (!Double.isNaN(correlation)) {
                    boolean flipalleles = wp.getFlipSNPAlleles()[d];
                    if (flipalleles) {
                        zscore = -zscore;
                        correlation = -correlation;
                    }
                    nrDatasetsPassingQC++;

                    double weight = Descriptives.getSqrt(numSamples);

                    zSum += (zscore * weight);
                    zSumAbsolute += (Math.abs(zscore) * weight);
                    nrTotalSamples += numSamples;
                    if (determinebeta) {
                        if (flipalleles) {
                            betasum += (-dsResults.beta[d][p] * numSamples);
                        } else {
                            betasum += (dsResults.beta[d][p] * numSamples);
                        }
                    }
                }
            }

            if (nrDatasetsPassingQC > 0 && nrTotalSamples > 0) {
                testsPerformed++;
                hasResults = true;
                double sqrtSample = Descriptives.getSqrt(nrTotalSamples);
                double zScore = zSum / sqrtSample;
                double zScoreAbsolute = zSumAbsolute / sqrtSample;
                double pValueOverall = Descriptives.convertZscoreToPvalue(zScore);
                double pValueOverallAbsolute = Descriptives.convertZscoreToPvalue(zScoreAbsolute);
                dsResults.pvalues[p] = pValueOverall;
                dsResults.pvaluesAbs[p] = pValueOverallAbsolute;
                dsResults.finalZScore[p] = zScore;
                dsResults.finalZScoreAbsolute[p] = zScoreAbsolute;
                // determine assessed allele....
                if (determinebeta) {
                    betasum /= nrTotalSamples;
                    double metase = 1 / Math.sqrt(nrTotalSamples);
                    dsResults.finalBeta[p] = betasum;
                    dsResults.finalBetaSe[p] = metase;
                }
            } else {
                dsResults.pvalues[p] = Double.NaN;
                dsResults.pvaluesAbs[p] = Double.NaN;
                dsResults.finalZScore[p] = Double.NaN;
                dsResults.finalZScoreAbsolute[p] = Double.NaN;
            }
            // calculate the weighted Z-score


        }

        wp.setHasResults(hasResults);

        wp.setResult(dsResults);

    }

    private void ploteQTL(WorkPackage wp, int p) {
        m_eQTLPlotter.draw(wp, p);
    }
    private int tmpbuffersize = 4096;

    private byte[] deflate(byte[] input) {
        Deflater d = new Deflater(6);
        d.setInput(input);
        d.finish();
        byte[] tmpbytebuffer = new byte[tmpbuffersize];
        int compressedDataLength = tmpbuffersize;
        int compressedsize = 0;
        byte[] finaldata = new byte[input.length + 1024];

        int start = 0;
        while (compressedDataLength == tmpbuffersize) {
            compressedDataLength = d.deflate(tmpbytebuffer);
            // out.write(bytebuffer, 0, compressedDataLength);

            System.arraycopy(tmpbytebuffer, 0, finaldata, start, compressedDataLength);
            start += compressedDataLength;
            compressedsize += compressedDataLength;
        }

        byte[] returndata = new byte[compressedsize];

        System.arraycopy(finaldata, 0, returndata, 0, compressedsize);
        return returndata;
    }

    private void deflateResults(WorkPackage currentWP) {
        Result r = currentWP.results;
//	double[][] datasetZscores = r.zscores;
        byte[][] inflatedZScores = new byte[r.zscores.length][0];
        if (r != null) {
            int[] numSamples = null;
            try {
                numSamples = r.numSamples;
            } catch (NullPointerException e) {
                System.out.println("ERROR: null result?");
            }

            double[][] zscores = r.zscores;
            int wpId = r.wpid;

            int[] probes = currentWP.getProbes();
            SNP[] snps = currentWP.getSnps();
            int numDatasets = zscores.length;

            if (zscores[0].length == 0) {
                System.out.println("Warning: Z-score list is empty!");
            }

            for (int d = 0; d < numDatasets; d++) {
                ByteBuffer buff = null;
                int nrBytesRequired = m_numProbes * 4;
                buff = ByteBuffer.allocate(nrBytesRequired);

                if (cisOnly) {
                    HashMap<Integer, Integer> availableProbes = new HashMap<Integer, Integer>();
                    int loc = 0;
                    for (Integer p : probes) {
                        availableProbes.put(p, loc);     // translate position if probe Id to probe[] (and zscore list) position
                        loc++;
                    }

                    for (int p = 0; p < m_numProbes; p++) {
                        Integer probeLoc = availableProbes.get(p);
                        if (probeLoc != null && !Double.isNaN(zscores[d][probeLoc])) {
                            if (currentWP.getFlipSNPAlleles()[d]) {
                                // zscorelist[p] = (float) -zscores[d][probeLoc];
                                buff.putFloat(p * 4, (float) -zscores[d][probeLoc]);
                            } else {
                                // zscorelist[p] = (float) zscores[d][probeLoc];
                                buff.putFloat(p * 4, (float) zscores[d][probeLoc]);
                            }
                        } else {
                            // zscorelist[p] = Float.NaN;
                            buff.putFloat(p * 4, Float.NaN);
                        }
                    }
                } else {
                    for (int i = 0; i < m_numProbes; i++) {
                        if (!Double.isNaN(zscores[d][i])) {
                            if (currentWP.getFlipSNPAlleles()[d]) {
                                // zscorelist[i] = (float) -zscores[d][i];
                                buff.putFloat(i * 4, (float) -zscores[d][i]);

                            } else {
                                // zscorelist[i] = (float) zscores[d][i];
                                buff.putFloat(i * 4, (float) zscores[d][i]);
                            }
                        } else {
                            // zscorelist[i] = Float.NaN;
                            buff.putFloat(i * 4, Float.NaN);
                        }
                    }
                }

                if (numSamples[d] != 0 && snps[d] != null) {
                    inflatedZScores[d] = deflate(buff.array());
                } else {
                    inflatedZScores[d] = null;
                }
            }
        }
        r.deflatedZScores = inflatedZScores;
    }
}
