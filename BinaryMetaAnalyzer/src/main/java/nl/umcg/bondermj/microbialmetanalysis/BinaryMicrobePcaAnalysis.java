/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.bondermj.microbialmetanalysis;

import nl.umcg.westrah.binarymetaanalyzer.*;
import gnu.trove.map.hash.TObjectIntHashMap;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Locale;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;

/**
 *
 * @author Marc Jan & Harm-Jan
 */
public class BinaryMicrobePcaAnalysis {

    public static void main(String[] args) {

        String settingsFile = null;
        String textToReplace = null;
        String replaceTextWith = null;

        if (args.length == 3) {
            settingsFile = args[0];
            textToReplace = args[1];
            replaceTextWith = args[2];

        } else if (args.length == 1) {
            settingsFile = args[0];

        } else {
            System.out.println("Usage: settings.xml replacetext replacetextwith");
            System.exit(-1);
        }

        BinaryMicrobePcaAnalysis meta = new BinaryMicrobePcaAnalysis(settingsFile, textToReplace, replaceTextWith);
        System.exit(0);

    }

    private MetaQTL4TraitAnnotation probeAnnotation;
    private BinaryMetaAnalysisDataset[] datasets = new BinaryMetaAnalysisDataset[0];
    private int[][] snpIndex;
    private String[] snpList;
    private final BinaryMetaAnalysisSettings settings;
    private String[] snpChr;
    private int[] snpPositions;
    private Integer[][] probeIndex;

    private QTL[] finalEQTLs;
    private boolean bufferHasOverFlown;
    private double maxSavedPvalue = -Double.MAX_VALUE;
    private int locationToStoreResult;

    TObjectIntHashMap<MetaQTL4MetaTrait> traitMap = new TObjectIntHashMap<MetaQTL4MetaTrait>();
    MetaQTL4MetaTrait[] traitList = null;

    public BinaryMicrobePcaAnalysis(String settingsFile, String textToReplace, String replaceTextWith) {
        // initialize settings
        settings = new BinaryMetaAnalysisSettings();
        settings.parse(settingsFile, textToReplace, replaceTextWith);
        int maxResults = settings.getFinalEQTLBufferMaxLength();
        int tmpbuffersize = (maxResults / 10);

        if (tmpbuffersize == 0) {
            tmpbuffersize = 10;
        } else if (tmpbuffersize > 250000) {
            tmpbuffersize = 250000;
        }

        try {
            run((maxResults + tmpbuffersize));
        } catch (IOException ex) {
            Logger.getLogger(BinaryMicrobePcaAnalysis.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    public void run(int bufferSize) throws IOException {

        String outdir = settings.getOutput();
        outdir = Gpio.formatAsDirectory(outdir);
        Gpio.createDir(outdir);
        // load probe annotation and index
        // this particular probe annotation can take multiple probes for a single location into account.
        System.out.println("Loading probe annotation from: " + settings.getProbetranslationfile());
        loadProbeAnnotation();

        for (int permutation = 0; permutation < settings.getNrPermutations() + 1; permutation++) {
            clearResultsBuffer();
            maxSavedPvalue = -Double.MAX_VALUE;
            // create dataset objects
            System.out.println("Running permutation " + permutation);
            datasets = new BinaryMetaAnalysisDataset[settings.getDatasetlocations().size()];

            System.out.println("Loading datasets");
            for (int d = 0; d < datasets.length; d++) {
                datasets[d] = new BinaryMetaAnalysisDataset(settings.getDatasetlocations().get(d), settings.getDatasetnames().get(d), settings.getDatasetPrefix().get(d), permutation, settings.getDatasetannotations().get(d), probeAnnotation);
            }
            System.out.println("Loaded " + datasets.length + " datasets");

            // create meta-analysis SNP index. have to recreate this every permutation, 
            // since the order of SNPs is generated at random.
            System.out.println("Creating SNP index");
            createSNPIndex();
            System.out.println("Total of " + snpIndex.length + " SNPs");
            System.out.println("Creating probe index");
            createProbeIndex();
            System.out.println("Total of " + probeIndex.length + " probes");

            if (snpChr == null) {
                System.out.println("Loading SNP annotation from " + settings.getSNPAnnotationFile());
                loadSNPAnnotation();
            }

            // run analysis
            System.out.println("Cis-analysis: " + settings.getAnalysisType());
            
            System.out.println("Cis-window: " + settings.getCisdistance());
            System.out.println("Trans-window: " + settings.getTransdistance());

            System.out.println("Starting meta-analysis");
            ProgressBar pb = new ProgressBar(snpList.length);

            snpLoop:
            for (int snp = 0; snp < snpList.length; snp++) {

                boolean printed = false;
                int[] sampleSizes = new int[datasets.length];
                int totalSampleSize = 0;
                Boolean[] flipZScores = new Boolean[datasets.length];
                String alleles = null;
                String alleleAssessed = null;

                // determine whether to flip the alleles for a certain dataset
                for (int d = 0; d < datasets.length; d++) {

                    int datasetSNPId = snpIndex[snp][d];
                    if (datasetSNPId == -9 && settings.isConfineSNPs()) {
                        continue snpLoop;
                    }
                    if (datasetSNPId != -9) {
                        sampleSizes[d] = datasets[d].getSampleSize(datasetSNPId);
                        totalSampleSize += sampleSizes[d];

                        if (alleles == null) {
                            flipZScores[d] = false;
                            alleles = datasets[d].getAlleles(datasetSNPId);
                            alleleAssessed = datasets[d].getAlleleAssessed(datasetSNPId);
                        } else {
                            String alleles2 = datasets[d].getAlleles(datasetSNPId);
                            String alleleAssessed2 = datasets[d].getAlleleAssessed(datasetSNPId);
                            flipZScores[d] = BaseAnnot.flipalleles(alleles, alleleAssessed, alleles2, alleleAssessed2);
                        }
                    }
                }
                Correlation.correlationToZScore(totalSampleSize + 10);

                Set<MetaQTL4MetaTrait> cisProbes = null;
                if (!settings.getAnalysisType().equals(BinaryMetaAnalysisSettings.Analysis.CISTRANS)) {
                    // do not test the cis probes
                    cisProbes = probeAnnotation.getMetatraits().getTraitInWindow(snpChr[snp], snpPositions[snp], settings.getTransdistance());
                }

                // iterate over the probe index
                float[][] finalZScores = new float[probeIndex.length][datasets.length];

                for (int d = 0; d < datasets.length; d++) {
                    if (datasets[d].getIsCisDataset()) {
                        System.err.println("ERROR: cannot run trans analysis on a cis dataset: " + settings.getDatasetlocations().get(d));
                        System.exit(-1);
                    }

                    if (flipZScores[d] == null) {
                        for (float[] zScore : finalZScores) {
                            zScore[d] = Float.NaN;
                        }
                    } else {
                        int datasetSNPId = snpIndex[snp][d];
                        float[] datasetZScores = datasets[d].getZScores(datasetSNPId);

// probeIndex[t.getMetaTraitId()][d] = p;
                        for (int p = 0; p < traitList.length; p++) {
                            MetaQTL4MetaTrait t = traitList[p];
                            if (cisProbes != null && cisProbes.contains(t)) {
                                finalZScores[p][d] = Float.NaN;
                            } else {
                                Integer datasetProbeId = probeIndex[p][d];
                                if (datasetProbeId != null) {
                                    finalZScores[p][d] = datasetZScores[datasetProbeId];
                                    if (flipZScores[d]) {
                                        finalZScores[p][d] *= -1;
                                    }
                                } else {
                                    finalZScores[p][d] = Float.NaN;
                                }
                            }
                        }
                    }
                }
                // meta-analyze!
                if (permutation > 0) {
                    double summedRsquare = 0;
                    for (int probe = 0; probe < traitList.length; probe++) {
                        double metaAnalysisZ = ZScores.getWeightedZ(finalZScores[probe], sampleSizes);
                        double tScore = ZScores.zScoreToCorrelation(metaAnalysisZ, totalSampleSize);
                        summedRsquare += tScore * tScore;
                    }
                    double newMetaZ = Correlation.convertCorrelationToZScore(totalSampleSize, Math.sqrt(summedRsquare));
                    double newMetaAnalysisP = Descriptives.convertZscoreToPvalue(newMetaZ);

                    // create output object
                    if (!Double.isNaN(newMetaAnalysisP)) {
                        MetaQTL4MetaTrait t = new MetaQTL4MetaTrait(21, "Microbe_Components", "-", -1, -1, "", traitList[0].getPlatformIds());
                        QTL q = new QTL(newMetaAnalysisP, t, snp, BaseAnnot.toByte(alleleAssessed), newMetaZ, BaseAnnot.toByteArray(alleles), new float[finalZScores[0].length], sampleSizes); // sort buffer if needed.
                        addEQTL(q);
                    } else {
                        System.out.println("Error in procedure.");
                    }
                } else {
                    double summedRsquare = 0;
                    float[] summedPerDataSet = new float[finalZScores[0].length];
                    for (int probe = 0; probe < traitList.length; probe++) {
                        double metaAnalysisZ = ZScores.getWeightedZ(finalZScores[probe], sampleSizes);
                        for (int i = 0; i < finalZScores[probe].length; i++) {
                            double tScore = ZScores.zScoreToCorrelation(finalZScores[probe][i], sampleSizes[i]);
                            summedPerDataSet[i] += tScore * tScore;
                        }
                        double tScore = ZScores.zScoreToCorrelation(metaAnalysisZ, totalSampleSize);
                        summedRsquare += tScore * tScore;
                    }

                    for (int i = 0; i < summedPerDataSet.length; i++) {
                        summedPerDataSet[i] = (float) Correlation.convertCorrelationToZScore(sampleSizes[i], Math.sqrt(summedPerDataSet[i]));
                    }
                    double newMetaZ = Correlation.convertCorrelationToZScore(totalSampleSize, Math.sqrt(summedRsquare));
                    double newMetaAnalysisP = Descriptives.convertZscoreToPvalue(newMetaZ);

                    // create output object
                    if (!Double.isNaN(newMetaAnalysisP)) {
                        MetaQTL4MetaTrait t = new MetaQTL4MetaTrait(21, "Microbe_Components", "-", -1, -1, "", traitList[0].getPlatformIds());
                        QTL q = new QTL(newMetaAnalysisP, t, snp, BaseAnnot.toByte(alleleAssessed), newMetaZ, BaseAnnot.toByteArray(alleles), summedPerDataSet, sampleSizes); // sort buffer if needed.
                        addEQTL(q);
                    } else {
                        System.out.println("Error in procedure.");
                    }
                }
                pb.iterate();
            }
            pb.close();

            for (BinaryMetaAnalysisDataset dataset : datasets) {
                dataset.close();
            }
            writeBuffer(outdir, permutation);
        }

        /*
         TODO:
         - ZSCORE RETRIEVAL
         - Plotting of z-scores
         - validation
         - multithreadalize
         */
    }

    private void createSNPIndex() throws IOException {

        HashSet<String> confineToTheseSNPs = null;
        if (settings.getSNPSelection() != null) {
            System.out.println("Selecting SNPs from file: " + settings.getSNPSelection());
            confineToTheseSNPs = new HashSet<String>();
            TextFile tf = new TextFile(settings.getSNPSelection(), TextFile.R);
            confineToTheseSNPs.addAll(tf.readAsArrayList());
            tf.close();

            System.out.println(confineToTheseSNPs.size() + " SNPs loaded.");
        }

        // create a list of all available SNPs
        HashSet<String> allSNPs = new HashSet<String>();
        for (BinaryMetaAnalysisDataset dataset : datasets) {
            String[] snps = dataset.getSNPs();
            for (String snp : snps) {
                if (confineToTheseSNPs == null || confineToTheseSNPs.contains(snp)) {
                    allSNPs.add(snp);
                }
            }

        }

        // create a temporary map that maps each SNP to a meta-analysis position
        int ctr = 0;
        TObjectIntHashMap<String> snpMap = new TObjectIntHashMap<String>(allSNPs.size(), 0.85f, -9);
        snpList = new String[allSNPs.size()];
        for (String s : allSNPs) {
            snpMap.put(s, ctr);
            snpList[ctr] = s;
            ctr++;
        }

        // fill index
        snpIndex = new int[allSNPs.size()][datasets.length];
        for (int d = 0; d < datasets.length; d++) {
            for (int s = 0; s < allSNPs.size(); s++) {
                snpIndex[s][d] = -9;
            }
        }
        for (int d = 0; d < datasets.length; d++) {
            String[] snps = datasets[d].getSNPs();
            for (int s = 0; s < snps.length; s++) {
                String snp = snps[s];
                int id = snpMap.get(snp);
                if (id != -9) {
                    snpIndex[id][d] = s;
                }
            }
        }
    }

    private void loadProbeAnnotation() throws IOException {

        HashSet<String> platforms = new HashSet<String>();
        platforms.addAll(settings.getDatasetannotations());
        probeAnnotation = new MetaQTL4TraitAnnotation(new File(settings.getProbetranslationfile()), platforms);
        traitList = new MetaQTL4MetaTrait[probeAnnotation.getMetatraits().size()];

        int q = 0;
        for (MetaQTL4MetaTrait t : probeAnnotation.getMetatraits()) {
            traitList[q] = t;
            traitMap.put(t, q);
            q++;
        }

    }

    private void loadSNPAnnotation() throws IOException {

        snpChr = new String[snpList.length];
        snpPositions = new int[snpList.length];
        for (int s = 0; s < snpList.length; s++) {
            snpChr[s] = "-10".intern();
            snpPositions[s] = -10;
        }

        TObjectIntHashMap<String> snpMap = new TObjectIntHashMap<String>(snpList.length);
        for (int s = 0; s < snpList.length; s++) {
            snpMap.put(snpList[s], s);
        }
        TextFile tf = new TextFile(settings.getSNPAnnotationFile(), TextFile.R);

        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String snp = elems[2];
            if (snpMap.contains(snp)) {
                int id = snpMap.get(snp);
                snpChr[id] = new String(elems[0].getBytes("UTF-8")).intern();
                snpPositions[id] = Integer.parseInt(elems[1]);
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

    }

    // index the probes 
    private void createProbeIndex() {
        probeIndex = new Integer[traitList.length][datasets.length];
        for (int d = 0; d < datasets.length; d++) {
            String[] probes = datasets[d].getProbeList();
            int platformId = probeAnnotation.getPlatformId(settings.getDatasetannotations().get(d));
            for (int p = 0; p < probes.length; p++) {
                MetaQTL4MetaTrait t = probeAnnotation.getTraitForPlatformId(platformId, probes[p]);
                int index = traitMap.get(t);
                probeIndex[index][d] = p;
            }
        }
    }

    private void addEQTL(QTL q) {

        double pval = q.getPvalue();
        if (bufferHasOverFlown) {
            if (pval <= maxSavedPvalue) {

                finalEQTLs[locationToStoreResult] = q;
                locationToStoreResult++;

                if (locationToStoreResult == finalEQTLs.length) {
                    Arrays.sort(finalEQTLs);
                    locationToStoreResult = settings.getFinalEQTLBufferMaxLength();
                    maxSavedPvalue = finalEQTLs[(settings.getFinalEQTLBufferMaxLength() - 1)].getPvalue();
                }
            }

        } else {
            if (pval > maxSavedPvalue) {
                maxSavedPvalue = pval;
            }

            finalEQTLs[locationToStoreResult] = q;
            locationToStoreResult++;

            if (locationToStoreResult == settings.getFinalEQTLBufferMaxLength()) {
                bufferHasOverFlown = true;
            }
        }
    }

    private void writeBuffer(String outdir, int permutation) throws IOException {

        // sort the finalbuffer for a last time
        if (locationToStoreResult != 0) {
            Arrays.sort(finalEQTLs, 0, locationToStoreResult);
        }

        String outfilename = outdir + "eQTLs.txt.gz";
        if (permutation > 0) {
            outfilename = outdir + "PermutedEQTLsPermutationRound" + permutation + ".txt.gz";
        }

        System.out.println("Writing output: " + outfilename);

        TextFile output = new TextFile(outfilename, TextFile.W);
        String header = "PValue\t"
                + "SNPName\t"
                + "SNPChr\t"
                + "SNPChrPos\t"
                + "ProbeName\t"
                + "ProbeChr\t"
                + "ProbeCenterChrPos\t"
                + "CisTrans\t"
                + "SNPType\t"
                + "AlleleAssessed\t"
                + "OverallZScore\t"
                + "DatasetsWhereSNPProbePairIsAvailableAndPassesQC\t"
                + "DatasetsZScores\t"
                + "DatasetsNrSamples\t"
                + "IncludedDatasetsMeanProbeExpression\t"
                + "IncludedDatasetsProbeExpressionVariance\t"
                + "HGNCName\t"
                + "IncludedDatasetsCorrelationCoefficient\t"
                + "Meta-Beta (SE)\t"
                + "Beta (SE)\t"
                + "FoldChange";

        output.writeln(header);
// PValue	SNPName	SNPChr	SNPChrPos	ProbeName	ProbeChr	ProbeCenterChrPos	CisTrans	SNPType	AlleleAssessed	OverallZScore	DatasetsWhereSNPProbePairIsAvailableAndPassesQC	DatasetsZScores	DatasetsNrSamples	IncludedDatasetsMeanProbeExpression	IncludedDatasetsProbeExpressionVariance	HGNCName	IncludedDatasetsCorrelationCoefficient	Meta-Beta (SE)	Beta (SE)	FoldChange	FDR

        DecimalFormat format = new DecimalFormat("###.#######", new DecimalFormatSymbols(Locale.US));
        DecimalFormat smallFormat = new DecimalFormat("0.#####E0", new DecimalFormatSymbols(Locale.US));
        for (int i = 0; i < settings.getFinalEQTLBufferMaxLength(); i++) {
            QTL q = finalEQTLs[i];
            if (q != null) {
                StringBuilder sb = new StringBuilder();
                if (q.getPvalue() < 1E-4) {
                    sb.append(smallFormat.format(q.getPvalue()));
                } else {
                    sb.append(format.format(q.getPvalue()));
                }

                sb.append("\t");
                int snpId = q.getSNPId();
                sb.append(snpList[snpId]);
                sb.append("\t");
                sb.append(snpChr[snpId]);
                sb.append("\t");
                sb.append(snpPositions[snpId]);
                sb.append("\t");
                MetaQTL4MetaTrait t = q.getMetaTrait();
                sb.append(t.getMetaTraitName());
                sb.append("\t");
                sb.append(t.getChr());
                sb.append("\t");
                sb.append(t.getChrMidpoint());
                sb.append("\t");
                int dist = Math.abs(t.getChrMidpoint() - snpPositions[snpId]);
                boolean sameChr = t.getChr().equals(snpChr[snpId]);
                if (sameChr && dist < settings.getCisdistance()) {
                    sb.append("Cis");
                } else if (sameChr && dist < settings.getTransdistance()) {
                    sb.append("Greyzone");
                } else {
                    sb.append("Trans");
                }

                sb.append("\t");
                sb.append(q.getAlleles());
                sb.append("\t");
                sb.append(q.getAlleleAssessed());
                sb.append("\t");
                sb.append(format.format(q.getZscore()));

                float[] datasetZScores = q.getDatasetZScores();
                String[] dsBuilder = new String[datasets.length];
                String[] dsNBuilder = new String[datasets.length];
                for (int d = 0; d < datasetZScores.length; d++) {
                    if (!Float.isNaN(datasetZScores[d])) {
                        dsBuilder[d] = settings.getDatasetnames().get(d);
                        dsNBuilder[d] = "" + q.getDatasetSampleSizes()[d];
                    } else {
                        dsBuilder[d] = "-";
                        dsNBuilder[d] = "-";
                    }
                }

                sb.append("\t");
                sb.append(Strings.concat(dsBuilder, Strings.semicolon));

                sb.append("\t");
                sb.append(Strings.concat(datasetZScores, format, Strings.semicolon));

                sb.append("\t");
                sb.append(Strings.concat(dsNBuilder, Strings.semicolon));
                sb.append("\t-\t-\t");

                sb.append(t.getAnnotation());
                sb.append("\t-\t-\t-\t-");

                output.writeln(sb.toString());
            }
        }

        output.close();

        System.out.println(
                "Done.");
    }

    private void clearResultsBuffer() {
        Arrays.fill(finalEQTLs, null);
        bufferHasOverFlown = false;
        locationToStoreResult = 0;
        maxSavedPvalue = -Double.MAX_VALUE;
    }
}
