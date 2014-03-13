package eqtlmappingpipeline.binarymeta.meta;

import eqtlmappingpipeline.binarymeta.meta.graphics.ZScorePlot;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.zip.DataFormatException;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.bin.BinaryResultDataset;
import umcg.genetica.math.stats.Descriptives;

public class ZScoreComparison extends MetaAnalyze {

    @Override
    public void analyze() throws IOException, DataFormatException, Exception {
        System.out.println("");
        System.out.println("Starting analysis!");

        String[] datasets = new String[m_settings.getDatasetnames().size()];
        for (int i = 0; i < m_settings.getDatasetnames().size(); i++) {
            datasets[i] = m_settings.getDatasetnames().get(i);
        }

        if (!m_settings.getOutput().endsWith("/")) {
            m_settings.setOutput(m_settings.getOutput() + "/MetaAnalysis/");
        }

        if (!Gpio.exists(m_settings.getOutput())) {
            Gpio.createDir(m_settings.getOutput());
        }

        String[] locations = new String[m_settings.getDatasetnames().size()];
        for (int i = 0; i < locations.length; i++) {
            locations[i] = m_settings.getDatasetlocations().get(i);
        }

        int permstart = 0;
        int permstop = m_settings.getNrPermutations() + 1;

        if (m_settings.getRunonlypermutation() > -1) {
            permstart = m_settings.getRunonlypermutation();
            permstop = m_settings.getRunonlypermutation() + 1;
        }

        for (int perm = permstart; perm < permstop; perm++) {
            ds = new BinaryResultDataset[m_settings.getDatasetlocations().size()];
            plotZScores(perm, locations, datasets);
        }

    }

    private void plotZScores(int perm, String[] locations, String[] datasets) throws IOException, Exception {
        pvaluedistribution = null;
        eQTLBuffer = null;
        finalEQTLBuffer = null;
        nrInFinalBuffer = 0;

        uniqueProbes = new HashSet<String>();
        uniqueSNPs = new HashSet<String>();

        int numDatasets = ds.length;
        probes = new ArrayList<String>();

        snps = new ArrayList<String>();
        snpChr = new ArrayList<Byte>();
        snpChrPos = new ArrayList<Integer>();

        nrTotalSamples = 0;
        String[] probeName = probeTranslation.getProbes();
        probes.addAll(Arrays.asList(probeName));

        initdatasets(locations, perm, -1);

        String zsName = null;
        if (m_settings.isMakezscoreplot()) {
            zs = new ZScorePlot();
            String[] datasets2 = new String[datasets.length + 1];
            System.arraycopy(datasets, 0, datasets2, 0, datasets.length);
            datasets2[datasets2.length - 1] = "Meta-Analysis";
            if (perm > 0) {
                zsName = m_settings.getOutput() + "ZScoreComparison-PermutationRound" + perm;
            } else {
                zsName = m_settings.getOutput() + "ZScoreComparison";
            }
            zs.init(numDatasets + 1, datasets2, true, zsName);
        }

        Descriptives.lookupSqrt(nrTotalSamples);
        pvaluedistribution = new int[m_settings.getNrOfBins()];

        eQTLBuffer = new EQTL[1];
        finalEQTLBuffer = new EQTL[0];

        pvaluethreshold = Double.MAX_VALUE;

        zsumPerSNP = new double[snps.size()];
        zsumSNPsNumberOfProbes = new int[snps.size()];
        zsumPerProbe = new double[probes.size()];
        zsumProbesNumberOfSNPs = new int[probes.size()];

        System.out.println("Performing the meta-analysis now: ");
//	System.out.println(snps.size() + "\t unique SNPs present in at least " + m_settings.snpDatasetPresenceThreshold + " datasets");
//	System.out.println(probes.size() + "\t unique Probespresent in at least " + m_settings.probeDatasetPresenceThreshold + " datasets");
        System.out.println(nrTotalSamples + "\t total samples");

//	if (m_settings.isMakezscoretable()) {
        if (perm == 0) {
            zscoretable = new TextFile(m_settings.getOutput() + "metazscoretable.txt.gz", TextFile.W, (10 * 1048576));
        } else {
            zscoretable = new TextFile(m_settings.getOutput() + "metazscoretable-Permutation" + perm + ".txt.gz", TextFile.W, (10 * 1048576));
        }
        StringBuilder zscoreout = new StringBuilder();
        zscoreout.append("SNP\tAlleleCoding\tAssessedAllele");
        for (int i = 0; i < probes.size(); i++) {
            zscoreout.append("\t").append(probes.get(i));
        }
        zscoretable.writeln(zscoreout.toString());

//	}

        /// init calculation pool,

        int nrProcs = Runtime.getRuntime().availableProcessors();
        if (m_settings.getNrThresds() > 0) {
            if (m_settings.getNrThresds() < nrProcs) {
                m_settings.setNrThresds(m_settings.getNrThresds());
            }
        }
        MetaAnalysisCalculationThread[] calcPool = new MetaAnalysisCalculationThread[nrProcs];
        LinkedBlockingQueue<MetaAnalysisWorkPackage> loaderQueue = new LinkedBlockingQueue<MetaAnalysisWorkPackage>(nrProcs);
        MetaAnalysisLoaderThread loaderThread = new MetaAnalysisLoaderThread(loaderQueue, snpTranslation, snps, ds);
        loaderThread.setName("Loader");
        loaderThread.start();

        PValueThreshold p = new PValueThreshold();
//	LinkedBlockingQueue<MetaAnalysisWorkPackage> resultQueue = new LinkedBlockingQueue<MetaAnalysisWorkPackage>(nrProcs);
//	MetaAnalysisResultThread resultThread = new MetaAnalysisResultThread(resultQueue, m_settings, datasets, perm, zscoretable, p);
//	resultThread.setName("Result");
//	resultThread.start();

        for (int i = 0; i < nrProcs; i++) {
            calcPool[i] = new MetaAnalysisPlotThread(loaderQueue, null, snps, probes, snpChr, snpChrPos, ds, snpTranslation, probeTranslationLookupTable, probeTranslation, m_settings, zs, p);
            calcPool[i].setName("MetaCalc-" + i);
            calcPool[i].start();
        }

        // kill the threads
        try {
            loaderThread.join();
            MetaAnalysisWorkPackage poison = new MetaAnalysisWorkPackage(0, 0);
            poison.poisonTheWell();

            for (int threadNum = 0; threadNum < calcPool.length; threadNum++) {
                try {
                    loaderQueue.put(poison);
                } catch (InterruptedException ex) {
                    ex.printStackTrace();
                }
            }
            for (int threadNum = 0; threadNum < calcPool.length; threadNum++) {
                calcPool[threadNum].join();
            }

//	    resultQueue.put(poison);
//	    resultThread.join();

        } catch (InterruptedException e) {
            System.err.println("Exception: Thread main interrupted.");
        }

        if (m_settings.isMakezscoretable()) {
//	    if (perm == 0) {
            zscoretable.close();
//	    }
        }
        if (zs != null) {
            zs.write(zsName);
        }
    }
}
