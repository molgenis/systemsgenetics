/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl3;

import eqtlmappingpipeline.metaqtl3.containers.EQTL;
import eqtlmappingpipeline.metaqtl3.containers.Result;
import eqtlmappingpipeline.metaqtl3.containers.WorkPackage;
import java.io.IOException;
import java.util.Arrays;
import java.util.concurrent.LinkedBlockingQueue;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDataset;
import umcg.genetica.io.trityper.eQTLTextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;

/**
 *
 * @author harmjan
 */
public class ResultProcessorThread extends Thread {

    
//    private BinaryResultProbeSummary[] m_dsProbeSummary = null;
//    private BinaryResultSNPSummary[] m_dsSNPSummary = null;
//    private BinaryGZipFloatMatrix[] m_dsZScoreMatrix = null;
//    private double m_pvaluethreshold = 2;
//    private int m_maxNrResults = 150000;
//    private final int m_totalNumberOfProbes;
//    private Result[] m_BinaryBuffer;
//    private final double m_pvaluePlotThreshold;
//    private final static char m_tab = '\t';
//    private int nrProcessed;
//    private EQTL[] tmpEQTLBuffer;
//    private int m_eQTLBufferCounter = 0;
//    private int m_result_counter = 0;
//    private int nrSet;
//    private int nrInFinalBuffer = 0;
//    private TextFile[] zScoreBinaryFile;
//    private TextFile zScoreMetaAnalysisFile; 
    
    
    private boolean m_createBinaryFiles = false;
    private TriTyperGeneticalGenomicsDataset[] m_gg = null;
    private boolean m_cisOnly;
    private Integer[][] m_probeTranslation;
    private int m_midpointprobedist;
    private final String m_outputdir;
    private final boolean m_permuting;
    private final int m_permutationround;
    private final boolean m_createTEXTFiles;
    private final String[] m_probeList;
    private final LinkedBlockingQueue<WorkPackage> m_queue;
    private final WorkPackage[] m_availableWorkPackages;
    private long nrTestsPerformed = 0;
    private EQTL[] finalEQTLs;
    private double maxSavedPvalue = Double.NaN;
    private int locationToStoreResult=0;
    private boolean bufferHasOverFlown=false;
    public int totalcounter = 0;
    private int m_maxResults = 0;
    private int m_numdatasets = 0;
    public double highestP = Double.MAX_VALUE;
    private int nrSNPsTested = 0;
    private final boolean m_useAbsoluteZScore;
    private BinaryFile[] zScoreBinaryFile;
    private BinaryFile zScoreMetaAnalysisFile;
    private TextFile zScoreMetaAnalysisRowNamesFile;
    private TextFile[] zScoreRowNamesFile;

    public ResultProcessorThread(int nrThreads, LinkedBlockingQueue<WorkPackage> queue, boolean chargeOutput,
            TriTyperGeneticalGenomicsDataset[] gg, MetaQTL3Settings settings, Integer[][] pprobeTranslation,
            boolean permuting, int round, String[] snplist, String[] probelist, WorkPackage[] allPackages) {
        m_availableWorkPackages = allPackages;
        m_createBinaryFiles = settings.createBinaryOutputFiles;
        m_createTEXTFiles = settings.createTEXTOutputFiles;
        m_useAbsoluteZScore = settings.useAbsoluteZScorePValue;
        m_queue = queue;
        m_outputdir = settings.outputReportsDir;
        m_permuting = permuting;
        m_permutationround = round;
        m_probeTranslation = pprobeTranslation;
        m_gg = gg;
        m_midpointprobedist = settings.ciseQTLAnalysMaxSNPProbeMidPointDistance;
        m_cisOnly = (settings.cisAnalysis && !settings.transAnalysis);

        m_probeList = probelist;
        m_maxResults = settings.maxNrMostSignificantEQTLs;

        int tmpbuffersize = m_maxResults / 10;
        if (tmpbuffersize == 0) {
            tmpbuffersize = 10;
        } 
        
//        else if (tmpbuffersize > 250000) {
//            tmpbuffersize = 250000;
//        }

//        m_totalNumberOfProbes = probelist.length;
//        m_pvaluePlotThreshold = settings.plotOutputPValueCutOff;
//        tmpEQTLBuffer = new EQTL[tmpbuffersize];
//        m_result_counter = 0;
        
        m_numdatasets = m_gg.length;

        finalEQTLs = new EQTL[(m_maxResults+tmpbuffersize)];
        nrSNPsTested = 0;
    }

    @Override
    public void run() {
//        nrProcessed = 0;
        try {
            if (m_createBinaryFiles) {
                zScoreBinaryFile = new BinaryFile[m_gg.length];
                zScoreRowNamesFile = new TextFile[m_gg.length];
                if (m_gg.length > 1) {
                    String metaAnalysisFileName = m_outputdir + "MetaAnalysis";
                    if (m_permuting) {
                        metaAnalysisFileName += "-PermutationRound-" + m_permutationround;
                    }
                    zScoreMetaAnalysisFile = new BinaryFile(metaAnalysisFileName + ".dat", BinaryFile.W);
                    zScoreMetaAnalysisRowNamesFile = new TextFile(metaAnalysisFileName + "-RowNames.txt.gz", TextFile.W);
                    zScoreMetaAnalysisRowNamesFile.writeln("SNP\tAlleles\tMinorAllele\tAlleleAssessed\tNrCalled");
                    TextFile tf = new TextFile(metaAnalysisFileName + "-ColNames.txt.gz", TextFile.W);
                    tf.writeList(Arrays.asList(m_probeList));
                    tf.close();
                }
                for (int d = 0; d < m_gg.length; d++) {
                    String fileName = m_outputdir + m_gg[d].getSettings().name;
                    if (m_permuting) {
                        fileName += "-PermutationRound-" + m_permutationround;
                    }
                    zScoreBinaryFile[d] = new BinaryFile(fileName + ".dat", BinaryFile.W);
                    TextFile tf = new TextFile(fileName + "-ColNames.txt.gz", TextFile.W);
                    tf.writeList(Arrays.asList(m_probeList));
                    tf.close();
                    zScoreRowNamesFile[d] = new TextFile(fileName + "-RowNames.txt.gz", TextFile.W);
                    zScoreRowNamesFile[d].writeln("SNP\tAlleles\tMinorAllele\tAlleleAssessed\tNrCalled\tMaf\tHWE\tCallRate");
                }
            }

            ProgressBar progressbar = new ProgressBar(m_availableWorkPackages.length);
            boolean poison = false;
            int counter = 0;
            while (!poison) {
                WorkPackage wp = m_queue.take();
                Result r = wp.results;
                if (wp.getHasResults()) {
                    nrSNPsTested++;
                }

                if (r.poison) {
                    poison = true;
                } else if (r.pvalues != null) {

                    nrTestsPerformed += wp.getNumTested();

                    
                    double[] pvalues = r.pvalues;
                    
                    //Is this working?
                    if (m_createBinaryFiles && !poison) {
                        writeBinaryResult(r);
                    }

                    if (m_createTEXTFiles && !poison) {
                        // classic textual output.

                        for (int p = 0; p < pvalues.length; p++) {
                            double pval = pvalues[p];

                            if (!Double.isNaN(pval) && pval <= highestP) {
                                double[][] corr = r.correlations;
                                double[] correlations = new double[corr.length];
                                double[] zscores = new double[corr.length];
                                int[] samples = new int[corr.length];

                                double[] fc = new double[corr.length];
                                double[] beta = new double[corr.length];
                                double[] betase = new double[corr.length];

                                for (int d = 0; d < correlations.length; d++) {
                                    if (Double.isNaN(corr[d][p])) {
                                        correlations[d] = Double.NaN;
                                        zscores[d] = Double.NaN;
                                        samples[d] = -9;
                                        fc[d] = Double.NaN;
                                        beta[d] = Double.NaN;
                                        betase[d] = Double.NaN;
                                    } else {
                                        correlations[d] = corr[d][p];
                                        if (m_useAbsoluteZScore) {
                                            zscores[d] = Math.abs(r.zscores[d][p]);
                                        } else {
                                            zscores[d] = r.zscores[d][p];
                                        }

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

            //Is this working?
            if (m_createBinaryFiles) {
                for (int d = 0; d < m_gg.length; d++) {
                    zScoreBinaryFile[d].close();
                    zScoreRowNamesFile[d].close();
                }
                if (m_gg.length > 1) {
                    zScoreMetaAnalysisFile.close();
                    zScoreMetaAnalysisRowNamesFile.close();
                }
            }

            if (m_createTEXTFiles) {
                if(locationToStoreResult>0){
                    if(!(bufferHasOverFlown&&locationToStoreResult==finalEQTLs.length)){
                        java.util.Arrays.sort(finalEQTLs);
                    }
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
           
            int wpId = r.wpid;
            WorkPackage currentWP = m_availableWorkPackages[wpId];
            double[][] zscores = r.zscores;
            if (m_cisOnly) {
                double[][] zscorestmp = new double[m_gg.length][m_probeList.length];
                int[] probes = currentWP.getProbes();
                for (int d = 0; d < zscorestmp.length; d++) {
                    for (int i = 0; i < zscorestmp[d].length; i++) {
                        zscorestmp[d][i] = Double.NaN;
                    }
                    for (int p = 0; p <probes.length; p++) {
                        int probeId = probes[p];
                        zscorestmp[d][probeId] = zscores[d][p]; 
                    }
                }
                zscores = zscorestmp;
            }

            if (zscores != null) {
                SNP[] snps = currentWP.getSnps();
                int numDatasets = zscores.length;
                double[] finalZscores = r.finalZScore;
                String snpoutput = null;
                if (m_gg.length > 1) {
                    int totalSampleNr = 0;
                    for (int d = 0; d < numDatasets; d++) {
                        if (snps[d] != null) {
                            String snpname = snps[d].getName();

                            byte[] alleles = snps[d].getAlleles();
                            byte minorAllele = snps[d].getMinorAllele();
                            byte alleleassessed = alleles[1];

                            if (currentWP.getFlipSNPAlleles()[d]) {
                                alleleassessed = alleles[0];
                            }
                            if (snpoutput == null) {
                                snpoutput = snpname + "\t" + BaseAnnot.getAllelesDescription(alleles) + "\t" + BaseAnnot.toString(minorAllele) + "\t" + BaseAnnot.toString(alleleassessed);
                            }
                            totalSampleNr += r.numSamples[d];
                        }
                    }
                    zScoreMetaAnalysisRowNamesFile.writeln(snpoutput + "\t" + totalSampleNr);

                    for (double z : finalZscores) {
                        zScoreMetaAnalysisFile.writeDouble(z);
                    }
                    
                }
                for (int d = 0; d < numDatasets; d++) {
                    double[] datasetZScores = zscores[d];
                    SNP datasetSNP = snps[d];
                    if (datasetSNP != null) {
                        BinaryFile outfile = zScoreBinaryFile[d];

                        String snpname = datasetSNP.getName();

                        byte[] alleles = datasetSNP.getAlleles();
                        byte minorAllele = datasetSNP.getMinorAllele();
                        byte alleleassessed = alleles[1];
                        double hwe = datasetSNP.getHWEP();
                        double cr = datasetSNP.getCR();
                        double maf = datasetSNP.getMAF();

                        if (currentWP.getFlipSNPAlleles()[d]) {
                            alleleassessed = alleles[0];
                        }
                        TextFile snpfile = zScoreRowNamesFile[d];
                        snpfile.writeln(snpname + "\t" + BaseAnnot.getAllelesDescription(alleles) + "\t" + BaseAnnot.toString(minorAllele) + "\t" + BaseAnnot.toString(alleleassessed) + "\t" + datasetSNP.getNrCalled() + "\t" + maf + "\t" + hwe + "\t" + cr);

                        for (double z : datasetZScores) {
                            outfile.writeDouble(z);
                        }
                        
                    }
                }
            }
        }
    }

    private void addEQTL(int pid, int sid, double pval, double zscore, double[] correlations, double[] zscores, int[] numSamples, byte[] alleles, byte assessedAllele, double[] fc, double[] beta, double[] betase, double finalbeta, double finalbetase) {

        if(bufferHasOverFlown){
            if(pval<maxSavedPvalue){
                EQTL e = new EQTL(m_numdatasets);
                e.pvalue = pval;
                e.pid = pid;
                e.sid = sid;
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
                finalEQTLs[locationToStoreResult] = e;
                locationToStoreResult++;
                totalcounter++;
            }
            if(locationToStoreResult==finalEQTLs.length){
                java.util.Arrays.sort(finalEQTLs);
                locationToStoreResult=m_maxResults;
                maxSavedPvalue = finalEQTLs[(m_maxResults-1)].pvalue;
            }
        } else {
            if(Double.isNaN(maxSavedPvalue)){
                maxSavedPvalue = pval;
            } else if(pval>maxSavedPvalue) {
                pval=maxSavedPvalue;
            }

            EQTL e = new EQTL(m_numdatasets);
            e.pvalue = pval;
            e.pid = pid;
            e.sid = sid;

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

            finalEQTLs[locationToStoreResult] = e;
            locationToStoreResult++;
            
            if(locationToStoreResult==m_maxResults){
                bufferHasOverFlown=true;
            }
        }
    }

//    private void mergeResults() {
//        EQTL[] toMerge = null;
//        if (m_eQTLBufferCounter < tmpEQTLBuffer.length) {
//            toMerge = new EQTL[m_eQTLBufferCounter];
//            System.arraycopy(tmpEQTLBuffer, 0, toMerge, 0, m_eQTLBufferCounter);
//        } else {
//            toMerge = tmpEQTLBuffer;
//        }
//
//        EQTL[] tmp = new EQTL[finalEQTLs.length + toMerge.length];
//        System.arraycopy(toMerge, 0, tmp, 0, toMerge.length);
//        System.arraycopy(finalEQTLs, 0, tmp, toMerge.length, finalEQTLs.length);
//
//        java.util.Arrays.sort(tmp);
//
//        nrInFinalBuffer += toMerge.length;
//        if (nrInFinalBuffer < m_maxResults) {
//            finalEQTLs = tmp;
//        } else {
//            finalEQTLs = new EQTL[m_maxResults];
//            System.arraycopy(tmp, 0, finalEQTLs, 0, m_maxResults);
//            for (int i = m_maxResults; i < tmp.length; i++) {
//                tmp[i].cleanUp();
//                tmp[i] = null;
//            }
//            nrInFinalBuffer = m_maxResults;
//            highestP = finalEQTLs[m_maxResults - 1].pvalue;
//        }
//    }

    private void writeTextResults() throws IOException {
        
        int nrOfEntriesToWrite = m_maxResults;
        if(locationToStoreResult<m_maxResults){
            nrOfEntriesToWrite = locationToStoreResult;
        }
        
        System.out.println("Writing " + nrOfEntriesToWrite + " results out of " + nrTestsPerformed + " tests performed. " + nrSNPsTested + " SNPs finally tested.");
        String fileName = m_outputdir + "eQTLs.txt.gz";
        if (m_permuting) {
            fileName = m_outputdir + "PermutedEQTLsPermutationRound" + m_permutationround + ".txt.gz";
            TextFile gz = new TextFile(fileName, TextFile.W);
            gz.writeln("PValue\tSNP\tProbe\tGene\tAlleles\tAlleleAssessed\tZScore");
            for (int i = 0; i < nrOfEntriesToWrite; i++) {
                String output = finalEQTLs[i].getDescription(m_availableWorkPackages, m_probeTranslation, m_gg, m_midpointprobedist);
                String[] realout = output.split("\t");
                String hugo = null;
                if (realout[eQTLTextFile.HUGO] == null) {
                    hugo = "-";
                } else {
                    hugo = realout[eQTLTextFile.HUGO];
                }
                String ln = realout[0] + "\t" + realout[1] + "\t" + realout[4] + "\t" + hugo + "\t" + realout[8] + "\t" + realout[9] + "\t" + realout[10];
                gz.writeln(ln);
            }
            gz.close();
        } else {
            eQTLTextFile et = new eQTLTextFile(fileName, eQTLTextFile.W);
            for (int i = 0; i < nrOfEntriesToWrite; i++) {
                String output = finalEQTLs[i].getDescription(m_availableWorkPackages, m_probeTranslation, m_gg, m_midpointprobedist);
                et.writeln(output);
            }
            et.close();
        }
    }
}
