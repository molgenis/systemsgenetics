/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl3;

import eqtlmappingpipeline.metaqtl3.containers.Settings;
import cern.colt.matrix.tint.IntMatrix2D;
import eqtlmappingpipeline.metaqtl3.containers.QTL;
import eqtlmappingpipeline.metaqtl3.containers.Result;
import eqtlmappingpipeline.metaqtl3.containers.WorkPackage;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.concurrent.LinkedBlockingQueue;
import javax.xml.bind.annotation.adapters.HexBinaryAdapter;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDataset;
import umcg.genetica.io.trityper.QTLTextFile;
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
//    private QTL[] tmpEQTLBuffer;
//    private int m_eQTLBufferCounter = 0;
//    private int m_result_counter = 0;
//    private int nrSet;
//    private int nrInFinalBuffer = 0;
//    private TextFile[] zScoreBinaryFile;
//    private TextFile zScoreMetaAnalysisFile; 
//    private int m_numdatasets = 0;
    long nrZ = 0;
    private boolean m_createBinaryFiles = false;
    private TriTyperGeneticalGenomicsDataset[] m_gg = null;
    private boolean m_cisOnly;
    private IntMatrix2D m_probeTranslation;
    private int m_midpointprobedist;
    private final String m_outputdir;
    private final boolean m_permuting;
    private final int m_permutationround;
    private final boolean m_createTEXTFiles;
    private final String[] m_probeList;
    private final LinkedBlockingQueue<WorkPackage> m_queue;
    private final WorkPackage[] m_availableWorkPackages;
    private long nrTestsPerformed = 0;
    private QTL[] finalEQTLs;
    private double maxSavedPvalue = -Double.MAX_VALUE;
    private int locationToStoreResult = 0;
    private boolean bufferHasOverFlown = false;
    private boolean sorted = false;
    private int m_maxResults = 0;
    public double highestP = Double.MAX_VALUE;
    private int nrSNPsTested = 0;
    private final boolean m_useAbsoluteZScore;
    private BinaryFile[] zScoreBinaryFile;
    private BinaryFile zScoreMetaAnalysisFile;
    private TextFile zScoreMetaAnalysisRowNamesFile;
    private TextFile[] zScoreRowNamesFile;

    public ResultProcessorThread(int nrThreads, LinkedBlockingQueue<WorkPackage> queue, boolean chargeOutput,
            TriTyperGeneticalGenomicsDataset[] gg, Settings settings, IntMatrix2D pprobeTranslation,
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

        int tmpbuffersize = (m_maxResults / 10);

        if (tmpbuffersize == 0) {
            tmpbuffersize = 10;
        } else if (tmpbuffersize > 250000) {
            tmpbuffersize = 250000;
        }

//        m_totalNumberOfProbes = probelist.length;
//        m_pvaluePlotThreshold = settings.plotOutputPValueCutOff;
//        tmpEQTLBuffer = new QTL[tmpbuffersize];
//        m_result_counter = 0;   
//        m_numdatasets = m_gg.length;
        finalEQTLs = new QTL[(m_maxResults + tmpbuffersize)];
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
                    // write magic number
                    if (m_cisOnly) {
                        zScoreMetaAnalysisFile.writeInt(1);
                    } else {
                        zScoreMetaAnalysisFile.writeInt(0);
                    }

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
                    // write magic number
                    if (m_cisOnly) {
                        zScoreBinaryFile[d].writeInt(1);
                    } else {
                        zScoreBinaryFile[d].writeInt(0);
                    }

                    TextFile tf = new TextFile(fileName + "-ColNames.txt.gz", TextFile.W);
                    tf.writeList(Arrays.asList(m_probeList));
                    tf.close();
                    zScoreRowNamesFile[d] = new TextFile(fileName + "-RowNames.txt.gz", TextFile.W);
                    zScoreRowNamesFile[d].writeln("SNP\tAlleles\tMinorAllele\tAlleleAssessed\tNrCalled\tMaf\tHWE\tCallRate");
                }
            }

            ProgressBar progressbar = new ProgressBar(m_availableWorkPackages.length);
            boolean poison = false;

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
            }

            progressbar.close();

            //Is this working?
            if (m_createBinaryFiles) {

                String fileName = "check";
                if (m_permuting) {
                    fileName += "-PermutationRound-" + m_permutationround;
                }
                fileName += ".md5";

                HexBinaryAdapter md5Parser = new HexBinaryAdapter();

                BufferedWriter md5writer = new BufferedWriter(new FileWriter(m_outputdir + fileName));

                for (int d = 0; d < m_gg.length; d++) {
                    zScoreBinaryFile[d].close();

                    fileName = m_gg[d].getSettings().name;
                    if (m_permuting) {
                        fileName += "-PermutationRound-" + m_permutationround;
                    }
                    fileName += ".dat";
                    md5writer.write(md5Parser.marshal(zScoreBinaryFile[d].getWrittenHash()) + "  " + fileName + '\n');

                    zScoreRowNamesFile[d].close();
                }
                if (m_gg.length > 1) {
                    zScoreMetaAnalysisFile.close();

                    fileName = "MetaAnalysis";
                    if (m_permuting) {
                        fileName += "-PermutationRound-" + m_permutationround;
                    }
                    fileName += ".dat";
                    md5writer.write(md5Parser.marshal(zScoreMetaAnalysisFile.getWrittenHash()) + "  " + fileName + '\n');

                    zScoreMetaAnalysisRowNamesFile.close();
                }

                md5writer.close();

            }

            if (m_createTEXTFiles) {
                if (!sorted) {
                    if (locationToStoreResult != 0) {

                        Arrays.parallelSort(finalEQTLs, 0, locationToStoreResult);
//                        SmoothSort.sort(finalEQTLs, 0, locationToStoreResult);
//                        inplaceArrayQuickSort.sort(finalEQTLs, 0, locationToStoreResult);

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

            if (zscores != null) {
                SNP[] snps = currentWP.getSnps();
                int numDatasets = zscores.length;
                double[] finalZscores = r.finalZScore;
                StringBuilder snpoutput = null;

                // if we're doing a meta-analysis, write the meta-analysis Z to a separate binaryFile
                if (m_gg.length > 1) {
                    int totalSampleNr = 0;
                    String snpname = null;
                    for (int d = 0; d < numDatasets; d++) {
                        if (snps[d] != null) {
                            snpname = snps[d].getName();

                            byte[] alleles = snps[d].getAlleles();
                            byte minorAllele = snps[d].getMinorAllele();
                            byte alleleassessed = alleles[1];

                            if (currentWP.getFlipSNPAlleles()[d]) {
                                alleleassessed = alleles[0];
                            }
                            if (snpoutput == null) {
                                snpoutput = new StringBuilder();
                                snpoutput.append(snpname);
                                snpoutput.append("\t");
                                snpoutput.append(BaseAnnot.getAllelesDescription(alleles));
                                snpoutput.append("\t");
                                snpoutput.append(BaseAnnot.toString(minorAllele));
                                snpoutput.append("\t");
                                snpoutput.append(BaseAnnot.toString(alleleassessed));
                            }
                            totalSampleNr += r.numSamples[d];
                        }
                    }

                    StringBuilder sb = null;
                    for (int p = 0; p < finalZscores.length; p++) {
                        float z = (float) finalZscores[p];
                        if (m_cisOnly) {
                            int[] probes = currentWP.getProbes();
                            int probeId = probes[p];
                            String probeName = m_probeList[probeId];
                            if (sb == null) {
                                sb = new StringBuilder();
                            } else {
                                sb.append("\t");
                            }
                            sb.append(probeName);

                            zScoreMetaAnalysisFile.writeFloat(z);
                        } else {
                            zScoreMetaAnalysisFile.writeFloat(z);
                        }
                    }

                    if (snpoutput != null) {
                        snpoutput.append("\t");
                        snpoutput.append(totalSampleNr);
                        snpoutput.append("\t-\t-\t-\t");
                        snpoutput.append(finalZscores.length);
                        snpoutput.append("\t");
                        if (sb != null) {
                            snpoutput.append(sb.toString());
                        } else {
                            snpoutput.append("-");
                        }
                        zScoreMetaAnalysisRowNamesFile.writeln(snpoutput.toString());
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
                        StringBuilder sb = null;
                        for (int p = 0; p < datasetZScores.length; p++) {
                            float z = (float) datasetZScores[p];
                            if (currentWP.getFlipSNPAlleles()[d]) {
                                z *= -1;
                            }
                            // System.out.println(p + "\t" + alleleassessed + "\t" + m_probeList[p] + "\t" + z + "\t" + currentWP.getFlipSNPAlleles()[d]);
                            if (m_cisOnly) {
                                // take into account that not all probes have been tested..
                                int[] probes = currentWP.getProbes();
                                int probeId = probes[p];
                                String probeName = m_probeList[probeId];
                                outfile.writeFloat(z);
                                if (sb == null) {
                                    sb = new StringBuilder();
                                } else {
                                    sb.append("\t");
                                }
                                sb.append(probeName);
                            } else {
                                outfile.writeFloat(z);
                            }
                        }

                        StringBuilder buffer = new StringBuilder();
                        buffer.append(snpname)
                                .append("\t")
                                .append(BaseAnnot.getAllelesDescription(alleles))
                                .append("\t")
                                .append(BaseAnnot.toString(minorAllele))
                                .append("\t")
                                .append(BaseAnnot.toString(alleleassessed))
                                .append("\t")
                                .append(datasetSNP.getNrCalled())
                                .append("\t")
                                .append(maf)
                                .append("\t")
                                .append(hwe)
                                .append("\t")
                                .append(cr)
                                .append("\t")
                                .append(datasetZScores.length)
                                .append("\t");
                        if (sb != null) {
                            buffer.append(sb.toString());
                        } else {
                            buffer.append("-");
                        }

                        snpfile.writeln(buffer.toString());

                    }
                }
            }
        }
    }

    private void addEQTL(int pid, int sid, double pval, double zscore, double[] correlations, double[] zscores, int[] numSamples, byte[] alleles, byte assessedAllele, double[] fc, double[] beta, double[] betase, double finalbeta, double finalbetase) {

        if (bufferHasOverFlown) {
            if (pval <= maxSavedPvalue) {

                sorted = false;

                finalEQTLs[locationToStoreResult] = new QTL(pval, pid, sid, assessedAllele, zscore, alleles, zscores, numSamples, correlations, fc, beta, betase, finalbeta, finalbetase);
                locationToStoreResult++;

                if (locationToStoreResult == finalEQTLs.length) {

                    Arrays.parallelSort(finalEQTLs);
//                    SmoothSort.sort(finalEQTLs);
//                    inplaceArrayQuickSort.sort(finalEQTLs);
                    sorted = true;
                    locationToStoreResult = m_maxResults;
                    maxSavedPvalue = finalEQTLs[(m_maxResults - 1)].getPvalue();
                }
            }

        } else {
            if (pval > maxSavedPvalue) {
                maxSavedPvalue = pval;
            }

            finalEQTLs[locationToStoreResult] = new QTL(pval, pid, sid, assessedAllele, zscore, alleles, zscores, numSamples, correlations, fc, beta, betase, finalbeta, finalbetase);
            locationToStoreResult++;

            if (locationToStoreResult == m_maxResults) {
                bufferHasOverFlown = true;
            }
        }
    }

    private void writeTextResults() throws IOException {

        int nrOfEntriesToWrite = m_maxResults;
        if (locationToStoreResult < m_maxResults) {
            nrOfEntriesToWrite = locationToStoreResult;
        }

        System.out.println("Writing " + nrOfEntriesToWrite + " results out of " + nrTestsPerformed + " tests performed. " + nrSNPsTested + " SNPs finally tested.");

        if (m_permuting) {
            TextFile gz = new TextFile((m_outputdir + "PermutedEQTLsPermutationRound" + m_permutationround + ".txt.gz"), TextFile.W);
            gz.writeln("PValue\tSNP\tProbe\tGene\tAlleles\tAlleleAssessed\tZScore");
            for (int i = 0; i < nrOfEntriesToWrite; i++) {
                gz.writeln(finalEQTLs[i].getPermutationDescription(m_availableWorkPackages, m_probeTranslation, m_gg, m_midpointprobedist));
            }
            gz.close();
        } else {
            QTLTextFile et = new QTLTextFile((m_outputdir + "eQTLs.txt.gz"), QTLTextFile.W);
            for (int i = 0; i < nrOfEntriesToWrite; i++) {
                et.writeln(finalEQTLs[i].getDescription(m_availableWorkPackages, m_probeTranslation, m_gg, m_midpointprobedist));
            }
            et.close();
        }
    }
}
