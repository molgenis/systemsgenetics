/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl4;

import java.util.Arrays;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.math.stats.Descriptives;

/**
 *
 * @author harmjan
 */
public class ResultProcessorThread extends Thread {

    private int[][] pvalueDistributions;
    private EQTL4[] finalEQTLs;
    private EQTL4[] tmpEQTLs;
    private int bufferSize = 50000;
    private EQTL4[][] probeLevelQTLs;
    private double pvalueThreshold = 2;
    private int nrProcessed;
    private int nrInBuffer = 0;
    private final int nrWorkPackagesToProcess;
    private final LinkedBlockingQueue<WorkPackage> queue;
    private final boolean performProbeLevelFDR;
    private int maxEQTLs;

    public ResultProcessorThread(int nrWorkPackagesToProcess, LinkedBlockingQueue<WorkPackage> queue, int[][] pvalueDistributions, EQTL4[] finalEQTLs, EQTL4[][] probeLevelQTLs, boolean performProbeLevelFDR, int maxEQTLs) {
        this.nrProcessed = 0;
        this.nrWorkPackagesToProcess = 0;
        this.queue = queue;
        this.finalEQTLs = finalEQTLs;
        this.probeLevelQTLs = probeLevelQTLs;
        this.performProbeLevelFDR = performProbeLevelFDR;
        this.tmpEQTLs = new EQTL4[bufferSize];
        this.maxEQTLs = maxEQTLs;
    }

    @Override
    public void run() {
        while (nrProcessed < nrWorkPackagesToProcess) {
            try {
                WorkPackage wp = queue.take();
                if (wp.getResults() != null) {
                    process(wp);
                }
                nrProcessed++;
            } catch (InterruptedException ex) {
                Logger.getLogger(ResultProcessorThread.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        if (nrInBuffer > 0) {
            mergeBuffers();
        }

        // now calculate FDR threshold
    }

    private void process(WorkPackage wp) {
        Integer[][] pvalueIndexes = wp.getResults();
        for (int perm = 0; perm < pvalueIndexes.length; perm++) {

            int[] pvalueDistribution = pvalueDistributions[0];
            if (perm > 0) {
                pvalueDistribution = pvalueDistributions[perm];
            }
            Integer[] vals = pvalueIndexes[perm];
            for (int probe = 0; probe < vals.length; probe++) {
                if (vals[probe] != null) {
                    Integer index = vals[probe];
                    pvalueDistribution[index]++;

                    double actualPvalue = Descriptives.m_zScoreToPValue[index];
                    int snp = wp.snp;
                    if (performProbeLevelFDR) {
                        EQTL4 eqtl = probeLevelQTLs[perm][probe];
                        if (eqtl == null) {
                            eqtl = new EQTL4(snp, probe, actualPvalue);
                        }

                        if (actualPvalue <= eqtl.getPvalue()) {
                            eqtl.setPvalue(actualPvalue);
                            eqtl.setSNP(snp);
                            eqtl.setProbe(probe);
                        }
                    }

                    // convert real data associations to eQTLs, and store them in the finalEQTL[]
                    if (perm == 0) {
                        if (actualPvalue < pvalueThreshold) {
                            tmpEQTLs[nrInBuffer] = new EQTL4(snp, probe, actualPvalue);
                            nrInBuffer++;

                            if (nrInBuffer == bufferSize) {
                                mergeBuffers();
                                nrInBuffer = 0;
                            }
                        }
                    }
                }
            }
        }
        wp.clearResults();
    }

    private void mergeBuffers() {
        boolean sort = true;
        if (finalEQTLs == null && nrInBuffer == bufferSize) {
            finalEQTLs = tmpEQTLs;
            pvalueThreshold = finalEQTLs[0].getPvalue();
        } else if (finalEQTLs == null && nrInBuffer < bufferSize) {
            finalEQTLs = new EQTL4[nrInBuffer];
            System.arraycopy(tmpEQTLs, 0, finalEQTLs, 0, nrInBuffer);
        } else if (finalEQTLs != null && finalEQTLs.length < maxEQTLs) {
            EQTL4[] tmp2EQTLs = new EQTL4[finalEQTLs.length + tmpEQTLs.length];
            System.arraycopy(finalEQTLs, 0, tmp2EQTLs, 0, finalEQTLs.length);
            System.arraycopy(tmp2EQTLs, 0, tmp2EQTLs, finalEQTLs.length, tmpEQTLs.length);
            finalEQTLs = tmp2EQTLs;
        } else {
            EQTL4[] tmp = new EQTL4[finalEQTLs.length + tmpEQTLs.length];
            System.arraycopy(tmpEQTLs, 0, tmp, 0, tmpEQTLs.length);
            System.arraycopy(finalEQTLs, 0, tmp, tmpEQTLs.length, finalEQTLs.length);
            Arrays.sort(tmp);
            finalEQTLs = new EQTL4[maxEQTLs];
            System.arraycopy(tmp, 0, finalEQTLs, 0, maxEQTLs);
            sort = false;
        }

        if (sort) {
            Arrays.sort(finalEQTLs);
        }
        tmpEQTLs = new EQTL4[bufferSize];
    }
}
