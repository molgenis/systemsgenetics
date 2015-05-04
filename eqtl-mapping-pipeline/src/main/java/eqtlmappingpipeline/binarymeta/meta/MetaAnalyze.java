///*
// * To change this template, choose Tools | Templates
// * and open the template in the editor.
// */
package eqtlmappingpipeline.binarymeta.meta;
//
//import eqtlmappingpipeline.gpio.binary.Dataset;

import eqtlmappingpipeline.binarymeta.meta.graphics.ZScorePlot;
import eqtlmappingpipeline.metaqtl3.FDR;
import eqtlmappingpipeline.metaqtl3.graphics.EQTLDotPlot;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.bin.BinaryResultDataset;
import umcg.genetica.io.trityper.bin.BinaryResultProbe;
import umcg.genetica.io.trityper.bin.BinaryResultSNP;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import umcg.genetica.math.stats.Descriptives;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.zip.DataFormatException;

///**
// *
// * @author harm-jan
// */
public class MetaAnalyze {

	protected static MetaSettings m_settings;
	protected BinaryResultDataset[] ds;
	protected ArrayList<String> probes;
	protected ArrayList<String> snps;
	protected Integer[][] snpTranslation;
	protected int[] pvaluedistribution;
	protected EQTL[] eQTLBuffer;
	protected EQTL[] finalEQTLBuffer;
	protected int nrInFinalBuffer = 0;
	protected double pvaluethreshold;
	protected ArrayList<Byte> snpChr;
	protected ArrayList<Integer> snpChrPos;
	protected ProbeTranslation probeTranslation;
	protected Integer[][] probeTranslationLookupTable;
	public static String header = "PValue\t"
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
			+ "IncludedDatasetsCorrelationCoefficient";
	protected double[] zsumPerSNP;
	protected int[] zsumSNPsNumberOfProbes;
	protected double[] zsumPerProbe;
	protected int[] zsumProbesNumberOfSNPs;
	protected ZScorePlot zs;
	protected TextFile zscoretable;
	protected HashSet<String> uniqueProbes;
	protected HashSet<String> uniqueSNPs;
	protected int nrTotalSamples;
	protected int numSNPs;
	protected int numProbes;
	private HashSet<String> probeListToAnalyze;

	public void init(String settingsFile, String texttoreplace, String replacetextwith) throws IOException {
		m_settings = new MetaSettings();
		m_settings.parse(settingsFile, texttoreplace, replacetextwith);
		probeTranslation = new ProbeTranslation();
		probeTranslation.load(m_settings.getProbetranslationfile());

	}

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
		m_settings.save();

		String[] locations = new String[m_settings.getDatasetnames().size()];
		for (int i = 0; i < locations.length; i++) {
			locations[i] = m_settings.getDatasetlocations().get(i);
		}

		int permstart = 0;
		int permstop = m_settings.getNrPermutations() + 1;

		if (m_settings.getRunonlypermutation() > -1) {
			permstart = m_settings.getRunonlypermutation();
			permstop = m_settings.getRunonlypermutation() + m_settings.getNrPermutations();
		}

		System.out.println(permstart + " - " + permstop);

		for (int perm = permstart; perm < permstop; perm++) {
			ds = new BinaryResultDataset[m_settings.getDatasetlocations().size()];
			runCalculationRound(perm, locations, datasets, -1);
		}

		if (m_settings.getRunonlypermutation() == -1) {

			if (m_settings.getNrPermutations() > 0) {
				FDR.calculateFDR(m_settings.getOutput(), m_settings.getNrPermutations(), m_settings.getFinalEQTLBufferMaxLength(), m_settings.getFdrthreshold(), true, null, null, FDR.FDRMethod.ALL, true);
				EQTLDotPlot edp = new EQTLDotPlot();
				edp.draw(m_settings.getOutput() + "/eQTLsFDR" + m_settings.getFdrthreshold() + ".txt", m_settings.getOutput() + "/DotPlot-FDR" + m_settings.getFdrthreshold() + ".pdf", EQTLDotPlot.Output.PDF); // "/eQTLsFDR" + fdrCutOff + ".txt", outputReportsDir + "/eQTLsFDR" + fdrCutOff + "DotPlot.png"
				edp = null;
			}
		}

	}

	protected void initdatasets(String[] locations, int perm, int dToUse) throws IOException {

		int numProbes = probeTranslation.getNumProbes();
		System.out.println(numProbes + " probes found in translation table. Now matching probes across datasets..");
		probeTranslationLookupTable = new Integer[ds.length][numProbes];
		HashSet<Integer> probesPresentInDatasets = new HashSet<Integer>();

//        m_settings.getSNPSelection();

		HashSet<String> selectedSNPs = null;

//        if (m_settings.getSNPSelection() != null) {
//            System.out.println("Selecting SNPs from: " + m_settings.getSNPSelection());
//            selectedSNPs = new HashSet<String>();
//            TextFile stf = new TextFile(m_settings.getSNPSelection(), TextFile.R);
//            selectedSNPs.addAll(stf.readAsArrayList());
//            stf.close();
//            System.out.println("Selected " + selectedSNPs.size() + " unique SNPs from file.");
//        }

		HashMap<String, HashSet<String>> selectedSNPProbePairs = null;
		if (m_settings.getSNPProbeSelection() != null) {
			System.out.println("Selecting SNP-probe pairs from: " + m_settings.getSNPProbeSelection());
			selectedSNPProbePairs = new HashMap<String, HashSet<String>>();
			selectedSNPs = new HashSet<String>();
			TextFile stf = new TextFile(m_settings.getSNPProbeSelection(), TextFile.R);
			int ctr = 0;
			String[] felems = stf.readLineElems(TextFile.tab);
			while (felems != null) {
				String snp = felems[0].intern();
				String probe = felems[1].intern();
				HashSet<String> probesForSNP = selectedSNPProbePairs.get(snp);
				if (probesForSNP == null) {
					probesForSNP = new HashSet<String>();
				}
				probesForSNP.add(probe.intern());
				selectedSNPs.add(snp.intern());
				selectedSNPProbePairs.put(snp.intern(), probesForSNP);
				ctr++;
				felems = stf.readLineElems(TextFile.tab);
			}

			stf.close();
			System.out.println("Selected " + ctr + " unique SNPs from file.");
		}

		HashSet<String> probesToInclude = null;

		if (m_settings.getProbeselection() != null) {
			TextFile tf = new TextFile(m_settings.getProbeselection(), TextFile.R);

			ArrayList<String> probesSelected = tf.readAsArrayList();

			probesToInclude = new HashSet<String>();
			probesToInclude.addAll(probesSelected);
			System.out.println(probesSelected.size() + " probes selected from file: " + m_settings.getProbeselection());
			tf.close();
		}

		for (int d = 0; d < ds.length; d++) {

			int probeAnnotationToUse = d;
			if (dToUse != -1) {
				probeAnnotationToUse = dToUse;
			}

			ds[d] = new BinaryResultDataset(locations[d], m_settings.getDatasetPrefix().get(probeAnnotationToUse), perm);
			BinaryResultProbe[] dsProbes = ds[d].getProbes();
			BinaryResultSNP[] dsSNPs = ds[d].getSnps();
			nrTotalSamples += ds[d].getMaxNrSamples();

			for (BinaryResultProbe p : dsProbes) {
				Integer newProbeId = probeTranslation.getProbeId(m_settings.getDatasetannotations().get(probeAnnotationToUse) + p.getName());
				if (newProbeId == null) {
					System.out.println(m_settings.getDatasetannotations().get(probeAnnotationToUse) + "\t" + p.getName() + " probe not present in annotationfile...?");
					System.exit(0);
				}
				if (probesToInclude == null || probesToInclude.contains("" + newProbeId)) {
					probesPresentInDatasets.add(newProbeId);
					probeTranslationLookupTable[d][newProbeId] = p.getId();
				} else {
					probeTranslationLookupTable[d][newProbeId] = null;
				}
			}

			for (BinaryResultSNP s : dsSNPs) {
				if (!uniqueSNPs.contains(s.getName().intern()) && (selectedSNPs == null || selectedSNPs.contains(s.getName().intern()))) {
					snps.add(s.getName().intern());
					snpChr.add(s.getChr());
					snpChrPos.add(s.getChrpos());
					uniqueSNPs.add(s.getName().intern());
				}
			}

//	    ds[d].clearProbeObjects();
		}

		TextFile probesPresentFile = new TextFile(m_settings.getOutput() + "ProbesPresentInAtLeastOneDataset.txt", TextFile.W);

		System.out.println(probesPresentInDatasets.size() + "\tunique probes present in all datasets.");
		Integer[] presentNrs = probesPresentInDatasets.toArray(new Integer[0]);
		for (Integer i : presentNrs) {
			probesPresentFile.writeln("" + i);
		}
		probesPresentFile.close();

		int selectedprobes = 0;


//        if (m_settings.getProbeselection() != null) {
//            TextFile tf = new TextFile(m_settings.getProbeselection(), TextFile.R);
//
//            String[] probesSelected = tf.readAsArray();
//
//            tf.close();
//
////	    for(int d=0; d<ds.length; d++){
////		Integer pid = probeTranslationLookupTable[d][6301];
////		if(pid!=null){
////		    System.out.println(6301+"\t"+pid+"\t"+ds[d].getProbes()[pid].getName());
////		} else {
////		    System.out.println(6301+"\t"+pid);
////		}
////
////	    }
////
////	    System.exit(0);
//            probeListToAnalyze = new HashSet<String>();
//            probeListToAnalyze.addAll(Arrays.asList(probesSelected));
//
//            System.out.println(probeListToAnalyze.size() + " unique probes selected for meta-analysis, from the file " + m_settings.getProbeselection());
//            int probePresenceCounter = 0;
//            for (int q = 0; q < probeTranslationLookupTable[0].length; q++) {
//                if (probeListToAnalyze != null) {
//                    if (!probeListToAnalyze.contains("" + q)) {
//
//                        for (int d = 0; d < ds.length; d++) {
//                            probeTranslationLookupTable[d][q] = null;
//                        }
//                    } else {
//                        probePresenceCounter++;
//                    }
//                }
//            }
//            System.out.println(probePresenceCounter + "\tprobes selected.");
//
//        } else {
		for (int q = 0; q < probeTranslationLookupTable[0].length; q++) {
			int probePresenceCounter = 0;
			if (probeListToAnalyze != null) {
				if (!probeListToAnalyze.contains("" + q)) {
					for (int d = 0; d < ds.length; d++) {
						probeTranslationLookupTable[d][q] = null;
					}
				}
			}

			for (int i = 0; i < ds.length; i++) {
				if (probeTranslationLookupTable[i][q] != null && ds[i].getMaxNrSamples() > m_settings.getProbeAndSNPPresenceFilterSampleThreshold()) {
					probePresenceCounter++;
				}
			}


			if (m_settings.getProbeDatasetPresenceThreshold() > 0 && probePresenceCounter < m_settings.getProbeDatasetPresenceThreshold()) {
				for (int d = 0; d < ds.length; d++) {
					probeTranslationLookupTable[d][q] = null;
				}
			} else if (probePresenceCounter > 0) {
				selectedprobes++;
			}
		}
		System.out.println("Selected " + selectedprobes + " probes from at least " + m_settings.getProbeDatasetPresenceThreshold() + " datasets of at least " + m_settings.getProbeAndSNPPresenceFilterSampleThreshold() + " samples.");
//        }
//	numProbes = uniqueProbes.size();
		numSNPs = uniqueSNPs.size();

		initSNPTranslation();
	}

	protected void initSNPTranslation() throws IOException {
		snpTranslation = new Integer[ds.length][numSNPs];

		for (int d = 0; d < ds.length; d++) {
			BinaryResultProbe[] dsProbes = ds[d].getProbes();
			BinaryResultSNP[] dsSNPs = ds[d].getSnps();

			for (int i = 0; i < snps.size(); i++) {
				BinaryResultSNP s = ds[d].getStringToSNP().get(snps.get(i));
				if (s != null) {
					snpTranslation[d][i] = s.getId();
				} else {
					snpTranslation[d][i] = null;
				}
			}
		}

		int selectedsnps = 0;

		HashSet<String> selectedSNPs = null;
		if (m_settings.getSNPSelection() != null) {
			System.out.println("Selecting SNPs from: " + m_settings.getSNPSelection());
			selectedSNPs = new HashSet<String>();
			TextFile stf = new TextFile(m_settings.getSNPSelection(), TextFile.R);
			selectedSNPs.addAll(stf.readAsArrayList());
			stf.close();
			System.out.println("Selected " + selectedSNPs.size() + " unique SNPs from file.");
		}

		TextFile selectedSNPFile = new TextFile(m_settings.getOutput() + "/SelectedSNPs.txt", TextFile.W);
		for (int s = 0; s < numSNPs; s++) {

			String snpName = snps.get(s);

			int snppresencecounter = 0;
			for (int d = 0; d < ds.length; d++) {
				if (snpTranslation[d][s] != null && ds[d].getMaxNrSamples() >= m_settings.getProbeAndSNPPresenceFilterSampleThreshold()) {
					snppresencecounter++;
				}
			}

			if (m_settings.getSnpDatasetPresenceThreshold() > 0 && snppresencecounter < m_settings.getSnpDatasetPresenceThreshold() || (selectedSNPs != null && !selectedSNPs.contains(snpName))) {
				for (int d = 0; d < ds.length; d++) {
					snpTranslation[d][s] = null;
				}
			} else if (snppresencecounter > 0) {
				selectedSNPFile.writeln(snps.get(s));
				selectedsnps++;
			}

		}

		selectedSNPFile.close();


		System.out.println("Selected " + selectedsnps + " snps from at least " + m_settings.getSnpDatasetPresenceThreshold() + " datasets of at least " + m_settings.getProbeAndSNPPresenceFilterSampleThreshold() + " samples.");
	}

	protected void runCalculationRound(int perm, String[] locations, String[] datasets, int dToUse) throws IOException, Exception {
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

		initdatasets(locations, perm, dToUse);

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

		eQTLBuffer = new EQTL[10000];
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

		if (m_settings.isMakezscoretable()) {
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

		}

		HashMap<String, HashSet<String>> selectedSNPProbePairs = null;
		if (m_settings.getSNPProbeSelection() != null) {
			System.out.println("Selecting SNP-probe pairs from: " + m_settings.getSNPProbeSelection());
			selectedSNPProbePairs = new HashMap<String, HashSet<String>>();

			TextFile stf = new TextFile(m_settings.getSNPProbeSelection(), TextFile.R);
			int ctr = 0;
			String[] felems = stf.readLineElems(TextFile.tab);
			while (felems != null) {
				String snp = felems[0];
				String probe = felems[1];
				HashSet<String> probesForSNP = selectedSNPProbePairs.get(snp);
				if (probesForSNP == null) {
					probesForSNP = new HashSet<String>();
				}
				probesForSNP.add(probe);
				selectedSNPProbePairs.put(snp, probesForSNP);
				ctr++;
				felems = stf.readLineElems(TextFile.tab);
			}

			stf.close();
			System.out.println("Selected " + ctr + " unique SNPs from file.");
		}

		/// init calculation pool,

		int nrProcs = Runtime.getRuntime().availableProcessors();
		if (m_settings.getNrThresds() > 0) {
			if (m_settings.getNrThresds() > nrProcs) {
				m_settings.setNrThresds(nrProcs);
			}
			nrProcs = m_settings.getNrThresds();
		}
		System.out.println("Using " + nrProcs + " threads :)");
		MetaAnalysisCalculationThread[] calcPool = new MetaAnalysisCalculationThread[nrProcs];
		LinkedBlockingQueue<MetaAnalysisWorkPackage> loaderQueue = new LinkedBlockingQueue<MetaAnalysisWorkPackage>(nrProcs);
		MetaAnalysisLoaderThread loaderThread = new MetaAnalysisLoaderThread(loaderQueue, snpTranslation, snps, ds);
		loaderThread.setName("Loader");
		loaderThread.start();

		PValueThreshold p = new PValueThreshold();
		LinkedBlockingQueue<MetaAnalysisWorkPackage> resultQueue = new LinkedBlockingQueue<MetaAnalysisWorkPackage>(nrProcs);
		MetaAnalysisResultThread resultThread = new MetaAnalysisResultThread(resultQueue, m_settings, datasets, perm, zscoretable, p, snps, selectedSNPProbePairs, probes);
		resultThread.setName("Result");
		resultThread.start();

		for (int i = 0; i < nrProcs; i++) {
			calcPool[i] = new MetaAnalysisCalculationThread(loaderQueue, resultQueue, snps, probes, snpChr, snpChrPos, ds, snpTranslation, probeTranslationLookupTable, probeTranslation, m_settings, zs, p);
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

			resultQueue.put(poison);
			resultThread.join();

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