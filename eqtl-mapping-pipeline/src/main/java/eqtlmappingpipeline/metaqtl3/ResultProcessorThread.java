/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl3;

import eqtlmappingpipeline.metaqtl3.containers.EQTL;
import eqtlmappingpipeline.metaqtl3.containers.Result;
import eqtlmappingpipeline.metaqtl3.containers.WorkPackage;
import eqtlmappingpipeline.metaqtl3.containers.eQTLResultContainer;
import java.io.IOException;
import java.util.concurrent.LinkedBlockingQueue;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDataset;
import umcg.genetica.io.trityper.bin.BinaryGZipFloatMatrix;
import umcg.genetica.io.trityper.bin.BinaryResultProbeSummary;
import umcg.genetica.io.trityper.bin.BinaryResultSNPSummary;
import umcg.genetica.io.trityper.eQTLTextFile;

/**
 *
 * @author harmjan
 */
public class ResultProcessorThread extends Thread {

    private boolean m_createBinaryFiles = false;
    private BinaryResultProbeSummary[] m_dsProbeSummary = null;
    private BinaryResultSNPSummary[] m_dsSNPSummary = null;
    private BinaryGZipFloatMatrix[] m_dsZScoreMatrix = null;
    private TriTyperGeneticalGenomicsDataset[] m_gg = null;
    private double m_pvaluethreshold = 2;
    //private int m_maxNrResults = 150000;
    private boolean m_cisOnly;
    private Integer[][] m_probeTranslation;
    private int m_midpointprobedist;
    private final String m_outputdir;
    private final boolean m_permuting;
    private final int m_permutationround;
    private final boolean m_createTEXTFiles;
    private final int m_totalNumberOfProbes;
    private final static char m_tab = '\t';
    private final String[] m_probeList;
    private Result[] m_BinaryBuffer;
    private final double m_pvaluePlotThreshold;
    private final LinkedBlockingQueue<WorkPackage> m_queue;
    private final WorkPackage[] m_availableWorkPackages;
    private int nrProcessed;
    private long nrTestsPerformed = 0;
    private EQTL[] tmpEQTLBuffer;
    private EQTL[] finalEQTLs;
    public int totalcounter = 0;
    private int m_eQTLBufferCounter = 0;
    private int m_result_counter = 0;
    private int m_maxResults = 0;
    private int m_numdatasets = 0;
    public double highestP = Double.MAX_VALUE;
    private int nrSet;
    private int nrInFinalBuffer = 0;

    public ResultProcessorThread(int nrThreads, LinkedBlockingQueue<WorkPackage> queue, boolean chargeOutput,
            TriTyperGeneticalGenomicsDataset[] gg, MetaQTL3Settings settings, Integer[][] pprobeTranslation,
            boolean permuting, int round, String[] snplist, String[] probelist, WorkPackage[] allPackages) {
        m_availableWorkPackages = allPackages;
        m_createBinaryFiles = settings.createBinaryOutputFiles;
        m_createTEXTFiles = settings.createTEXTOutputFiles;
        m_queue = queue;
        m_outputdir = settings.outputReportsDir;
        m_totalNumberOfProbes = probelist.length;
        m_permuting = permuting;
        m_permutationround = round;
        m_probeTranslation = pprobeTranslation;
        m_gg = gg;
        m_midpointprobedist = settings.ciseQTLAnalysMaxSNPProbeMidPointDistance;
        m_cisOnly = (settings.cisAnalysis && !settings.transAnalysis);

        m_probeList = probelist;
        m_maxResults = settings.maxNrMostSignificantEQTLs;
        m_pvaluePlotThreshold = settings.plotOutputPValueCutOff;

        int tmpbuffersize = m_maxResults / 4;
        if (tmpbuffersize == 0) {
            tmpbuffersize = 10;
        } else if (tmpbuffersize > 250000) {
            tmpbuffersize = 250000;
        }
        tmpEQTLBuffer = new EQTL[tmpbuffersize];

        m_result_counter = 0;
        m_numdatasets = m_gg.length;

        finalEQTLs = new EQTL[0];
    }

    @Override
    public void run() {
        nrProcessed = 0;
        try {
            if (m_createBinaryFiles) {
                m_dsProbeSummary = new BinaryResultProbeSummary[m_gg.length];
                m_dsSNPSummary = new BinaryResultSNPSummary[m_gg.length];
                m_dsZScoreMatrix = new BinaryGZipFloatMatrix[m_gg.length];

                for (int d = 0; d < m_gg.length; d++) {
                    String fileName = m_outputdir + m_gg[d].getSettings().name;
                    if (m_permuting) {
                        fileName += "-PermutationRound-" + m_permutationround;
                    } else {
                        BinaryResultProbeSummary dsProbeSummary = new BinaryResultProbeSummary(fileName, BinaryResultProbeSummary.W);
                        dsProbeSummary.write(m_probeList, m_probeTranslation, m_gg);
                        dsProbeSummary.close();
                    }

                    m_dsSNPSummary[d] = new BinaryResultSNPSummary(fileName, BinaryResultSNPSummary.W);
                    m_dsZScoreMatrix[d] = new BinaryGZipFloatMatrix(fileName, BinaryGZipFloatMatrix.W);
                }
            }

            ProgressBar progressbar = new ProgressBar(m_availableWorkPackages.length);
            boolean poison = false;
            int counter = 0;
            while (!poison) {
                WorkPackage wp = m_queue.take();

                Result r = wp.results;

                if (r.poison) {
                    poison = true;
                } else if (r.pvalues != null) {

                    nrTestsPerformed += wp.getNumTested();
                    
                    double[] pvalues = r.pvalues;
                    if (m_createBinaryFiles && !poison) {
                        writeBinaryResult(r);
                    }

                    if (m_createTEXTFiles && !poison) {
                        // classic textual output.

                        for (int p = 0; p < pvalues.length; p++) {
                            double pval = pvalues[p];
//                            if (!Double.isNaN(pval)) {
//                                nrTestsPerformed++;
//                            }
                            if (!Double.isNaN(pval) && pval <= highestP) {
                                double[][] corr = r.correlations;
                                Double[] correlations = new Double[corr.length];
                                Double[] zscores = new Double[corr.length];
                                Integer[] samples = new Integer[corr.length];

                                Double[] fc = new Double[corr.length];
                                Double[] beta = new Double[corr.length];
                                Double[] betase = new Double[corr.length];

                                for (int d = 0; d < correlations.length; d++) {
                                    if (Double.isNaN(corr[d][p])) {
                                        correlations[d] = null;
                                        zscores[d] = null;
                                        samples[d] = null;
                                        fc[d] = null;
                                        beta[d] = null;
                                        betase[d] = null;
                                    } else {
                                        correlations[d] = corr[d][p];
                                        zscores[d] = r.zscores[d][p];
                                        samples[d] = r.numSamples[d];
                                        fc[d] = r.fc[d][p];
                                        beta[d] = r.beta[d][p];
                                        betase[d] = r.se[d][p];
                                    }
                                }
//
                                byte allele = -1;
                                byte[] alleles = null;
                                SNP[] snps = wp.getSnps();
                                for (int d = 0; d < snps.length; d++) {
                                    if (snps[d] != null) {
                                        allele = snps[d].getMinorAllele();
                                        alleles = snps[d].getAlleles();
                                        break;
                                    }
                                }


                                if (alleles == null) {
                                    System.err.println("SNP has null alleles: ");
                                    for (int d = 0; d < snps.length; d++) {

                                        if (snps[d] != null) {

                                            allele = snps[d].getMinorAllele();
                                            System.err.println(allele);
                                            alleles = snps[d].getAlleles();
                                            System.err.println(alleles);
                                            break;
                                        }
                                    }
                                }


                                double Zfinal = r.finalZScore[p];
                                double finalbeta = r.finalBeta[p];
                                double finalbetase = r.finalBetaSe[p];
                                int pid;
                                if (m_cisOnly) {
                                    pid = wp.getProbes()[p];
                                } else {
                                    pid = p;
                                }

                                addEQTL(pid, wp.getId(), pval, Zfinal, correlations, zscores, samples, alleles, allele, fc, beta, betase, finalbeta, finalbetase);

                            }
                        }
                    }

                }

                if (wp.results != null) {
                    wp.clearResults();

                }

                progressbar.iterate();
                counter++;
            }

            progressbar.close();

            if (m_createBinaryFiles) {
                for (int d = 0; d < m_gg.length; d++) {
                    m_dsSNPSummary[d].close();
                    m_dsZScoreMatrix[d].flush();
                    m_dsZScoreMatrix[d].close();
                }
            }

            if (m_createTEXTFiles) {
                if (m_eQTLBufferCounter > 0) {
                    mergeResults();
                }
                writeTextResults();
            }

        } catch (IOException e1) {
            e1.printStackTrace();
        } catch (InterruptedException e2) {
            e2.printStackTrace();
        }
    }

    private void writeBinaryResult(Result r) throws IOException {
        if (r != null) {
            int[] numSamples = null;
            try {
                numSamples = r.numSamples;
            } catch (NullPointerException e) {
                System.out.println("ERROR: null result?");
            }

            double[][] zscores = r.zscores;
            int wpId = r.wpid;

            WorkPackage currentWP = m_availableWorkPackages[wpId];

            int[] probes = currentWP.getProbes();
            SNP[] snps = currentWP.getSnps();
            int numDatasets = zscores.length;

            if (zscores[0].length == 0) {
                System.out.println("Warning: Z-score list is empty!");
            }

            if (r.deflatedZScores == null) {
                System.err.println("Error: deflated binary data expected, but not detected.");
                System.exit(-1);
            }
            for (int d = 0; d < numDatasets; d++) {
                byte[] buff = r.deflatedZScores[d];
                if (buff != null && numSamples[d] != 0 && snps[d] != null) {
                    long index = m_dsZScoreMatrix[d].writeDeflated(buff);
                    String snpname = snps[d].getName();
                    byte snpchr = snps[d].getChr();
                    int snpchrpos = snps[d].getChrPos();
                    double hwe = snps[d].getHWEP();
                    double cr = snps[d].getCR();
                    double maf = snps[d].getMAF();
                    byte[] alleles = snps[d].getAlleles();
                    byte minorAllele = snps[d].getMinorAllele();
                    byte alleleassessed = alleles[1];

                    if (currentWP.getFlipSNPAlleles()[d]) {
                        alleleassessed = alleles[0];
                    }

                    m_dsSNPSummary[d].write(snpname, snpchr, snpchrpos, hwe, maf, cr, alleles, minorAllele, alleleassessed, numSamples[d], index);
                }
            }
        }
    }

    private void addEQTL(int pid, int sid, double pval, double zscore, Double[] correlations, Double[] zscores, Integer[] numSamples, byte[] alleles, byte assessedAllele, Double[] fc, Double[] beta, Double[] betase, double finalbeta, double finalbetase) {
        EQTL e = new EQTL(m_numdatasets);
        e.pvalue = pval;
        e.pid = pid;
        e.sid = sid;
        //if (!m_permuting) {
        e.alleleAssessed = assessedAllele;
        e.zscore = zscore;
        e.alleles = alleles;
        e.datasetZScores = zscores;
        e.datasetsSamples = numSamples;
        e.correlations = correlations;
        e.datasetfc = fc;
        e.datasetbeta = beta;
        e.datasetbetase = betase;
        e.finalbeta = finalbeta;
        e.finalbetase = finalbetase;
        //}
        tmpEQTLBuffer[m_eQTLBufferCounter] = e;

        m_eQTLBufferCounter++;
        totalcounter++;

        if (m_eQTLBufferCounter == tmpEQTLBuffer.length) {
            mergeResults();
            m_eQTLBufferCounter = 0;
        }
    }

    private void mergeResults() {
        EQTL[] toMerge = null;
        if (m_eQTLBufferCounter < tmpEQTLBuffer.length) {
            toMerge = new EQTL[m_eQTLBufferCounter];
            System.arraycopy(tmpEQTLBuffer, 0, toMerge, 0, m_eQTLBufferCounter);
        } else {
            toMerge = tmpEQTLBuffer;
        }


        EQTL[] tmp = new EQTL[finalEQTLs.length + toMerge.length];
        System.arraycopy(toMerge, 0, tmp, 0, toMerge.length);
        System.arraycopy(finalEQTLs, 0, tmp, toMerge.length, finalEQTLs.length);

        java.util.Arrays.sort(tmp);

        nrInFinalBuffer += toMerge.length;
        if (nrInFinalBuffer < m_maxResults) {
            finalEQTLs = tmp;
        } else {
            finalEQTLs = new EQTL[m_maxResults];
            System.arraycopy(tmp, 0, finalEQTLs, 0, m_maxResults);
            for (int i = m_maxResults; i < tmp.length; i++) {
                tmp[i].cleanUp();
                tmp[i] = null;
            }
            nrInFinalBuffer = m_maxResults;
            highestP = finalEQTLs[m_maxResults - 1].pvalue;
        }
    }

    private void writeTextResults() throws IOException {
        System.out.println("Writing " + finalEQTLs.length + " results out of " + nrTestsPerformed);
        String fileName = m_outputdir + "eQTLs.txt.gz";
        if (m_permuting) {
            fileName = m_outputdir + "PermutedEQTLsPermutationRound" + m_permutationround + ".txt.gz";
            TextFile gz = new TextFile(fileName, TextFile.W);
            gz.writeln("PValue\tSNP\tProbe\tGene");
            for (int i = 0; i < finalEQTLs.length; i++) {
                String output = finalEQTLs[i].getDescription(m_availableWorkPackages, m_probeTranslation, m_gg, m_midpointprobedist);
                String[] realout = output.split("\t");
                String hugo = null;
                if (realout[eQTLTextFile.HUGO] == null) {
                    hugo = "-";
                } else {
                    hugo = realout[eQTLTextFile.HUGO];
                }
                String ln = realout[0] + "\t" + realout[1] + "\t" + realout[4] + "\t" + hugo;
                gz.writeln(ln);
            }
            gz.close();
        } else {
            eQTLTextFile et = new eQTLTextFile(fileName, eQTLTextFile.W);
            for (int i = 0; i < finalEQTLs.length; i++) {
                String output = finalEQTLs[i].getDescription(m_availableWorkPackages, m_probeTranslation, m_gg, m_midpointprobedist);
                et.writeln(output);
            }
            et.close();
        }
    }
}
