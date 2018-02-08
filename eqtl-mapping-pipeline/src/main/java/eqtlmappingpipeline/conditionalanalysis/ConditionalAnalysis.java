/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.conditionalanalysis;

import eqtlmappingpipeline.metaqtl3.EQTLRegression;
import eqtlmappingpipeline.metaqtl3.FDR;
import eqtlmappingpipeline.metaqtl3.MetaQTL3;
import gnu.trove.set.hash.THashSet;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import umcg.genetica.console.ConsoleGUIElems;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDataset;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDatasetSettings;
import umcg.genetica.io.trityper.QTLTextFile;
import umcg.genetica.io.trityper.util.ChrAnnotation;
import umcg.genetica.io.trityper.util.DetermineLD;
import umcg.genetica.text.Strings;

/**
 * @author harmjan
 */
public class ConditionalAnalysis extends MetaQTL3 {
	
	private boolean skipinitialSNPMapping = false;
	private boolean skipAlleQTLMapping = false;
	
	public void setSkipInitialSNPMapping() {
		skipinitialSNPMapping = true;
	}
	
	public void run(String xmlSettingsFile, String texttoreplace, String texttoreplacewith,
					String ingt, String inexp, String inexpplatform, String inexpannot, String gte,
					String out, boolean cis, boolean trans, int perm, boolean textout, boolean binout, String snpfile, Integer threads) throws IOException, Exception {
		
		initialize(xmlSettingsFile, texttoreplace, texttoreplacewith, null, null, ingt, inexp, inexpplatform, inexpannot, gte, out, cis, trans, perm, textout, binout, snpfile, threads, null, null, null, true, true, null, null, null);
		double fdrthreshold = m_settings.fdrCutOff;
		
		m_settings.provideBetasAndStandardErrors = true;
		m_settings.provideFoldChangeData = true;
		String origOutputDir = m_settings.outputReportsDir;
		m_settings.outputReportsDir = origOutputDir + "SNP-Initial" + Gpio.getFileSeparator();
		m_settings.plotOutputDirectory = origOutputDir + "SNP-Initial" + Gpio.getFileSeparator() + "plots" + Gpio.getFileSeparator();
		Gpio.createDir(m_settings.plotOutputDirectory);
		Gpio.createDir(m_settings.outputReportsDir);
		
		String[] gzfiles = Gpio.getListOfFiles(origOutputDir, "gz");
		for (String filename : gzfiles) {
			File file = new File(filename);
			String f = file.getName();
			Gpio.moveFile(filename, m_settings.outputReportsDir + f);
		}
		
		System.out.println(m_settings.outputReportsDir + " is used for output");
		
		if (!skipinitialSNPMapping && !skipAlleQTLMapping) {
			// run the eQTL analysis over the SNPs --> determine strongest probes for each SNP
			mapEQTLs();
		}
		
		HashSet<String> tsSNPConfine = m_settings.tsSNPsConfine;
		
		// take the significant eQTL Probes
		String fdrSignificantFile = "eQTLsFDR" + fdrthreshold + "-ProbeLevel.txt";
		String fdrAllFile = "eQTLsFDR-ProbeLevel.txt.gz";
		
		FDR.FDRMethod fdrType = m_settings.fdrType;
		if (fdrType.equals(FDR.FDRMethod.GENELEVEL)) {
			fdrSignificantFile = "eQTLsFDR" + fdrthreshold + "-GeneLevel.txt";
			fdrAllFile = "eQTLsFDR-GeneLevel.txt.gz";
		} else if (fdrType.equals(FDR.FDRMethod.ALL)) {
			fdrSignificantFile = "eQTLsFDR" + fdrthreshold + "-ProbeLevel.txt";
			fdrAllFile = "eQTLsFDR-ProbeLevel.txt.gz";
		} else if (fdrType.equals(FDR.FDRMethod.FULL)) {
			fdrSignificantFile = "eQTLsFDR" + fdrthreshold + ".txt";
			fdrAllFile = "eQTLsFDR.txt.gz";
		}
		
		QTLTextFile etf = new QTLTextFile(origOutputDir + "/SNP-Initial/" + fdrSignificantFile, QTLTextFile.R);
		EQTL[] eQTLsSNPsUnconditional = etf.read();
		etf.close();
		
		System.out.println("Loaded: " + eQTLsSNPsUnconditional.length + " eQTL from " + origOutputDir + "/SNP-Initial/" + fdrSignificantFile);
		
		// run the eQTL analysis over the probes -- > get strongest effect per probe
		m_settings.outputReportsDir = origOutputDir + "/Probe-Initial/";
		m_settings.plotOutputDirectory = origOutputDir + "/Probe-Initial/plots/";
		Gpio.createDir(m_settings.outputReportsDir);
		Gpio.createDir(m_settings.plotOutputDirectory);
		m_settings.tsSNPsConfine = null;
		m_settings.tsProbesConfine = new THashSet<String>();
		
		HashMap<String, HashSet<String>> esnpProbePairs = new HashMap<String, HashSet<String>>();
		HashMap<String, HashSet<String>> gsnpProbePairs = new HashMap<String, HashSet<String>>();
		
		for (EQTL e : eQTLsSNPsUnconditional) {
			m_settings.tsProbesConfine.add(e.getProbe());
			HashSet<String> probes = gsnpProbePairs.get(e.getRsName());
			if (probes == null) {
				probes = new HashSet<String>();
			}
			probes.add(e.getProbe());
			gsnpProbePairs.put(e.getRsName(), probes);
		}
		for (TriTyperGeneticalGenomicsDatasetSettings s : m_settings.datasetSettings) {
			s.tsProbesConfine = m_settings.tsProbesConfine;
		}
		
		System.out.println(m_settings.tsProbesConfine.size() + "\tprobes have significant effects.");
		if (!skipAlleQTLMapping) {
			reinit(null);
			mapEQTLs();
		}
		
		etf = new QTLTextFile(origOutputDir + "/Probe-Initial/" + fdrAllFile, QTLTextFile.R);
		EQTL[] eQTLsProbesUnconditional = etf.read();
		etf.close();
		
		HashSet<String> eSNPsConfine = new HashSet<String>();
		
		HashMap<String, EQTL> probeeQTLMap = new HashMap<String, EQTL>();
		for (EQTL e : eQTLsProbesUnconditional) {
			if (!probeeQTLMap.containsKey(e.getProbe())) {
				probeeQTLMap.put(e.getProbe(), e);
				HashSet<String> probes = esnpProbePairs.get(e.getRsName());
				if (probes == null) {
					probes = new HashSet<String>();
				}
				eSNPsConfine.add(e.getRsName());
				probes.add(e.getProbe());
				esnpProbePairs.put(e.getRsName(), probes);
			}
		}
		
		//eQTLsProbesUnconditional
		// run conditional analysis: test only SNP/Probe pairs, condition on top probe effect
		m_settings.outputReportsDir = origOutputDir + "/SNPs-Conditional/";
		m_settings.plotOutputDirectory = origOutputDir + "/SNPs-Conditional/plots/";
		Gpio.createDir(m_settings.plotOutputDirectory);
		Gpio.createDir(m_settings.outputReportsDir);
		
		ArrayList<EQTL> eQTLsToRegress = new ArrayList<EQTL>();
		HashSet<String> probes = new HashSet<String>();
		for (EQTL e : eQTLsProbesUnconditional) {
			if (!probes.contains(e.getProbe())) {
				eQTLsToRegress.add(e);
				probes.add(e.getProbe());
			}
		}
		EQTL[] eQTLsToRegressArr = eQTLsToRegress.toArray(new EQTL[0]);

//        m_settings.regressOutEQTLEffectFileName = origOutputDir + "/Probe-Initial/eQTLProbesFDR" + fdrthreshold + "-ProbeLevel.txt";
		m_settings.tsSNPsConfine = tsSNPConfine;
		m_settings.tsSNPProbeCombinationsConfine = gsnpProbePairs;
		m_settings.performEQTLAnalysisOnSNPProbeCombinationSubset = true;
		if (!skipAlleQTLMapping) {
			reinit(eQTLsToRegressArr);
			mapEQTLs();
		}
		// run conditional analysis: test ony SNP/Probe pairs, condition on snp effect
		m_settings.outputReportsDir = origOutputDir + "/Probes-Conditional/";
		m_settings.plotOutputDirectory = origOutputDir + "/Probes-Conditional/plots/";
		Gpio.createDir(m_settings.plotOutputDirectory);
		Gpio.createDir(m_settings.outputReportsDir);
//        m_settings.regressOutEQTLEffectFileName = origOutputDir + "/SNP-Initial/eQTLProbesFDR" + fdrthreshold + "-ProbeLevel.txt";
		
		eQTLsToRegress = new ArrayList<EQTL>();
		probes = new HashSet<String>();
		for (EQTL e : eQTLsSNPsUnconditional) {
			if (!probes.contains(e.getProbe())) {
				eQTLsToRegress.add(e);
				probes.add(e.getProbe());
			}
		}
		eQTLsToRegressArr = eQTLsToRegress.toArray(new EQTL[0]);
		
		m_settings.tsSNPsConfine = eSNPsConfine;
		m_settings.tsSNPProbeCombinationsConfine = esnpProbePairs;
		m_settings.performEQTLAnalysisOnSNPProbeCombinationSubset = true;
		if (!skipAlleQTLMapping) {
			reinit(eQTLsToRegressArr);
			mapEQTLs();
		}
		
		System.out.println(ConsoleGUIElems.DOUBLELINE);
		System.out.println("Done with eQTL mappings.. Now summarizing results.");
		
		// now summarize the results in some clever way.
		HashMap<Pair<String, String>, Double> fdrSNPsConditional = readFDRFile(origOutputDir + "/SNPs-Conditional/" + fdrAllFile);
		HashMap<Pair<String, String>, Double> fdrProbesConditional = readFDRFile(origOutputDir + "/Probes-Conditional/" + fdrAllFile);
		HashMap<Pair<String, String>, EQTL> snpConditionaleQTLMap = readEQTLFile(origOutputDir + "/SNPs-Conditional/" + fdrAllFile);
		HashMap<Pair<String, String>, EQTL> probeConditionaleQTLMap = readEQTLFile(origOutputDir + "/Probes-Conditional/" + fdrAllFile);

//        HashMap<Pair<String, String>, Double> fdrSNPsConditional = readFDRFile(origOutputDir + "/SNPs-Conditional/eQTLsFDR.txt.gz");
//        HashMap<Pair<String, String>, Double> fdrSNPsConditional = readFDRFile(origOutputDir + "/SNPs-Conditional/eQTLsFDR.txt.gz");
		DetermineLD ld = new DetermineLD();
		
		TextFile outfile = new TextFile(origOutputDir + "/ConditionalAnalysis.txt", TextFile.W);
		outfile.writeln("Probe\tProbe Chr\tProbe ChrPos\tProbe HUGO\t"
				+ "gSNP\tgSNP Chr\tgSNP ChrPos\tgSNP alleles\tgSNP allele assessed\t"
				+ "eSNP\teSNP Chr\teSNP ChrPos\teSNP alleles\teSNP allele assessed\t"
				+ "LD (rsquared, dprime)\t"
				+ "eSNP N\tgSNP N\t"
				+ "gSNP-PVal\teSNP-PVal\tgSNP-PVal-Conditional\teSNP-PVal-Conditional\t"
				+ "gSNP-ZMeta\teSNP-ZMeta\tgSNP-ZMeta-Conditional\teSNP-ZMeta-Conditional\t"
				+ "gSNP-BMeta\teSNP-BMeta\tgSNP-BMeta-Conditional\teSNP-BMeta-Conditional\t"
				+ "gSNP-B\teSNP-B\tgSNP-B-Conditional\teSNP-B-Conditional\t"
				//                + "gSNP-FC\teSNP-FC\tgSNP-FC-Conditional\teSNP-FC-Conditional\t"
				+ "gSNP-FDR\teSNP-FDR\tgSNP-FDR-Conditional\teSNP-FDR-Conditional\t");
		
		SNPLoader[] loaders = new SNPLoader[m_gg.length];
		for (int d = 0; d < m_gg.length; d++) {
			loaders[d] = m_gg[d].getGenotypeData().createSNPLoader();
		}
		
		for (EQTL snpUnconditionalEQTL : eQTLsSNPsUnconditional) {
			EQTL probeunconditionaleqtl = probeeQTLMap.get(snpUnconditionalEQTL.getProbe());
			if (probeunconditionaleqtl == null) {
				System.err.println("ERROR: probe " + snpUnconditionalEQTL.getProbe() + " was not tested in probe-unconditional analysis.");
				System.exit(0);
			}
			
			Pair<String, String> gSNPPair = new Pair<String, String>(snpUnconditionalEQTL.getRsName(), snpUnconditionalEQTL.getProbe());
			Pair<String, String> eSNPPair = new Pair<String, String>(probeunconditionaleqtl.getRsName(), probeunconditionaleqtl.getProbe());
//	    System.out.println(e.getRsName() + "\t" + e.getProbe() + "\t" + probeunconditionaleqtl.getRsName());
			EQTL snpconditionaleqtl = snpConditionaleQTLMap.get(gSNPPair);
			EQTL probeconditionaleqtl = probeConditionaleQTLMap.get(eSNPPair);
			
			String gsnp = snpUnconditionalEQTL.getRsName();
			String esnp = probeunconditionaleqtl.getRsName();
			String ldstr = "";
			String genename = snpUnconditionalEQTL.getProbeHUGO();
			if (genename == null) {
				genename = "-";
			}
			String output = snpUnconditionalEQTL.getProbe() + "\t" + ChrAnnotation.parseByte(snpUnconditionalEQTL.getProbeChr()) + "\t" + snpUnconditionalEQTL.getProbeChrPos() + "\t" + genename
					+ "\t" + gsnp + "\t" + ChrAnnotation.parseByte(snpUnconditionalEQTL.getRsChr()) + "\t" + snpUnconditionalEQTL.getRsChrPos() + "\t" + snpUnconditionalEQTL.getAlleles() + "\t" + snpUnconditionalEQTL.getAlleleAssessed()
					+ "\t" + esnp + "\t" + ChrAnnotation.parseByte(probeunconditionaleqtl.getRsChr()) + "\t" + probeunconditionaleqtl.getRsChrPos() + "\t" + probeunconditionaleqtl.getAlleles() + "\t" + probeunconditionaleqtl.getAlleleAssessed() + "\t";
			
			if (gsnp.equals(esnp)) {
				ldstr = "1";
			} else {
				String[] ldStrArr = new String[m_gg.length];
				for (int ds = 0; ds < m_gg.length; ds++) {
					ldStrArr[ds] = "-";
					TriTyperGeneticalGenomicsDataset d = m_gg[ds];
					Integer gsnpid = d.getGenotypeData().getSnpToSNPId().get(gsnp);
					Integer esnpid = d.getGenotypeData().getSnpToSNPId().get(esnp);
					if (gsnpid != -9 && esnpid != -9) {
						SNP gsnpobj = d.getGenotypeData().getSNPObject(gsnpid);
						SNP esnpobj = d.getGenotypeData().getSNPObject(esnpid);
						if (gsnpobj != null && esnpobj != null) {
							loaders[ds].loadGenotypes(gsnpobj);
							loaders[ds].loadGenotypes(esnpobj);
							double r2 = ld.getRSquared(gsnpobj, esnpobj, d.getGenotypeData(), DetermineLD.RETURN_R_SQUARED, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);
							double dprime = ld.getRSquared(gsnpobj, esnpobj, d.getGenotypeData(), DetermineLD.RETURN_D_PRIME, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);
							ldStrArr[ds] = "" + r2 + ", " + dprime;
						}
						gsnpobj.clearGenotypes();
						esnpobj.clearGenotypes();
					}
				}
				ldstr = Strings.concat(ldStrArr, Strings.semicolon);
			}
			
			Integer[] samples = snpUnconditionalEQTL.getDatasetsSamples();
			int sampleSizeSNPUnconditional = 0;
			for (int i = 0; i < samples.length; i++) {
				if (samples[i] != null) {
					sampleSizeSNPUnconditional += samples[i];
				}
			}
			
			samples = probeunconditionaleqtl.getDatasetsSamples();
			int sampleSizeProbeUnconditional = 0;
			for (int i = 0; i < samples.length; i++) {
				if (samples[i] != null) {
					sampleSizeProbeUnconditional += samples[i];
				}
			}
			
			String eSNPSample = "" + sampleSizeSNPUnconditional;
			String gSNPSample = "" + sampleSizeProbeUnconditional;
			String probeunconditionalpval = "" + probeunconditionaleqtl.getPvalue();
			String snpconditionalpval = "-";
			String snpconditionalz = "-";
			String snpconditionalmetabeta = "-";
			String snpconditionalbeta = "-";
			String snpconditionalfc = "-";
			String snpconditionalfdr = "-";
			
			String probeconditionalpval = "-";
			String probeconditionalz = "-";
			String probeconditionalmetabeta = "-";
			String probeconditionalbeta = "-";
			String probeconditionalfc = "-";
			String probeconditionalfdr = "-";
			
			if (snpconditionaleqtl != null) {
				snpconditionalpval = "" + snpconditionaleqtl.getPvalue();
				snpconditionalmetabeta = snpconditionaleqtl.getMetaBeta();
				snpconditionalbeta = snpconditionaleqtl.getBeta();
				snpconditionalfc = snpconditionaleqtl.getFC();
				snpconditionalz = "" + snpconditionaleqtl.getZscore();
				snpconditionalfdr = "" + fdrSNPsConditional.get(gSNPPair);
			}
			
			if (probeconditionaleqtl != null) {
				probeconditionalpval = "" + probeconditionaleqtl.getPvalue();
				probeconditionalmetabeta = probeconditionaleqtl.getMetaBeta();
				probeconditionalbeta = probeconditionaleqtl.getBeta();
				probeconditionalfc = probeconditionaleqtl.getFC();
				probeconditionalz = "" + probeconditionaleqtl.getZscore();
				probeconditionalfdr = "" + fdrProbesConditional.get(eSNPPair);
			}
			
			output += ldstr
					+ "\t" + eSNPSample
					+ "\t" + gSNPSample
					+ "\t" + snpUnconditionalEQTL.getPvalue()
					+ "\t" + probeunconditionalpval
					+ "\t" + snpconditionalpval
					+ "\t" + probeconditionalpval
					+ "\t" + snpUnconditionalEQTL.getZscore()
					+ "\t" + probeunconditionaleqtl.getZscore()
					+ "\t" + snpconditionalz
					+ "\t" + probeconditionalz
					+ "\t" + snpUnconditionalEQTL.getMetaBeta()
					+ "\t" + probeunconditionaleqtl.getMetaBeta()
					+ "\t" + snpconditionalmetabeta
					+ "\t" + probeconditionalmetabeta
					+ "\t" + snpUnconditionalEQTL.getBeta()
					+ "\t" + probeunconditionaleqtl.getBeta()
					+ "\t" + snpconditionalbeta
					+ "\t" + probeconditionalbeta
					+ "\t" + snpUnconditionalEQTL.getFDR()
					+ "\t" + probeunconditionaleqtl.getFDR()
					+ "\t" + snpconditionalfdr
					+ "\t" + probeconditionalfdr;
			
			outfile.writeln(output);
			System.out.println(output);
			
		}
		for (int ds = 0; ds < m_gg.length; ds++) {
			loaders[ds].close();
			loaders[ds].close();
		}
		outfile.close();
	}
	
	private void reinit(EQTL[] regressthis) throws IOException, Exception {
		m_gg = null;
		
		int numDatasets = m_settings.datasetSettings.size();
		m_gg = new TriTyperGeneticalGenomicsDataset[numDatasets];
		numAvailableInds = 0;
		for (int i = 0; i < numDatasets; i++) {
			
			System.out.println("- Loading dataset: " + m_settings.datasetSettings.get(i).name + "");
			System.out.println(ConsoleGUIElems.LINE);
			m_gg[i] = new TriTyperGeneticalGenomicsDataset(m_settings.datasetSettings.get(i));
		}
		
		for (int i = 0; i < numDatasets; i++) {
			if (!m_settings.performParametricAnalysis) {
				
				m_gg[i].getExpressionData().rankAllExpressionData(m_settings.equalRankForTies);
			}
			m_gg[i].getExpressionData().calcAndSubtractMean();
			m_gg[i].getExpressionData().calcMeanAndVariance();
			numAvailableInds += m_gg[i].getExpressionToGenotypeIdArray().length;
			
		}
		
		if (regressthis != null && regressthis.length > 0) {
			EQTLRegression eqr = new EQTLRegression();
			eqr.regressOutEQTLEffects(regressthis, m_gg);
			// rerank
			for (int i = 0; i < numDatasets; i++) {
				if (!m_settings.performParametricAnalysis) {
					
					m_gg[i].getExpressionData().rankAllExpressionData(m_settings.equalRankForTies);
				}
				m_gg[i].getExpressionData().calcAndSubtractMean();
				m_gg[i].getExpressionData().calcMeanAndVariance();
				numAvailableInds += m_gg[i].getExpressionToGenotypeIdArray().length;
				
			}
		}
		
		
		System.out.println(ConsoleGUIElems.LINE);
		System.out.println("");
		
		// sort datasets on size to increase efficiency of random reads..
		// for some reason, it is faster to load the largest dataset first.
		Arrays.sort(m_gg, Collections.reverseOrder());
		
		System.out.println("Accumulating available data...");
		System.out.print(ConsoleGUIElems.LINE);
		
		createSNPList();
		createProbeList();
		
		// create WorkPackage objects
		determineSNPProbeCombinations();
		
		if (m_workPackages == null || m_workPackages.length == 0) {
			System.err.println("Error: No work detected");
			System.exit(0);
		}
		
		// determine number of threadss
		if (m_settings.nrThreads == null) {
			m_settings.nrThreads = Runtime.getRuntime().availableProcessors();
		} else {
			int numProcs = Runtime.getRuntime().availableProcessors();
			if (m_settings.nrThreads > numProcs || m_settings.nrThreads < 1) {
				m_settings.nrThreads = numProcs;
			}
		}
		
		if (m_workPackages.length < m_settings.nrThreads) {
			m_settings.nrThreads = m_workPackages.length;
		}
		printSummary();
	}
	
	public void runGetResultsForAllSNPs(String xmlSettingsFile, String texttoreplace, String texttoreplacewith,
										String ingt, String inexp, String inexpplatform, String inexpannot, String gte,
										String out, boolean cis, boolean trans, int perm, boolean textout, boolean binout, String snpfile, Integer threads) throws IOException, Exception {
		// initialize..
//        initialize(xmlSettingsFile, texttoreplace, texttoreplacewith, ingt, inexp, inexpplatform, inexpannot, gte, out, cis, trans, perm, textout, binout, snpfile, threads, null, null, null);
//
//        double fdrthreshold = m_settings.fdrCutOff;
//        m_settings.provideBetasAndStandardErrors = true;
//        m_settings.provideFoldChangeData = true;
//        String origOutputDir = m_settings.outputReportsDir;
//        m_settings.outputReportsDir = origOutputDir + "/SNP-Initial/";
//        m_settings.plotOutputDirectory = origOutputDir + "/SNP-Initial/plots/";
//        Gpio.createDir(m_settings.plotOutputDirectory);
//        Gpio.createDir(m_settings.outputReportsDir);
//
//        String[] gzfiles = Gpio.getListOfFiles(origOutputDir, "gz");
//        for (String filename : gzfiles) {
//            String f = filename.replace(origOutputDir, "");
//            Gpio.moveFile(filename, m_settings.outputReportsDir + f);
//        }
//
//
//        System.out.println(m_settings.outputReportsDir + " is used for output");
//
//        if (!skipinitialSNPMapping) {
//            // run the eQTL analysis over the SNPs --> determine strongest probes for each SNP
//            mapEQTLs();
//        }
//
//        HashSet<String> tsSNPConfine = m_settings.tsSNPsConfine;
//
//        // the result is a set of SNP probe combo's. This is a SNP-centric analysis, so we want to load the top probe per SNP, irrespective of the significance.
//        HashSet<String> probesToTest = new HashSet<String>();
//        TextFile tf = new TextFile(origOutputDir + "/SNP-Initial/eQTLs.txt", TextFile.R);
//        tf.readLine();
//        String[] elems = tf.readLineElems(TextFile.tab);
//        while (elems != null) {
//            probesToTest.add(elems[4]);
//            elems = tf.readLineElems(TextFile.tab);
//        }
//        tf.close();
//
//
//        // run the eQTL analysis over the probes -- > get strongest effect per probe
//        m_settings.outputReportsDir = origOutputDir + "/Probe-Initial/";
//        m_settings.plotOutputDirectory = origOutputDir + "/Probe-Initial/plots/";
//        Gpio.createDir(m_settings.outputReportsDir);
//        Gpio.createDir(m_settings.plotOutputDirectory);
//        m_settings.tsSNPsConfine = null;
//        m_settings.tsProbesConfine = probesToTest;
//        for (TriTyperGeneticalGenomicsDatasetSettings s : m_settings.datasetSettings) {
//            s.tsProbesConfine = m_settings.tsProbesConfine;
//        }
//        System.out.println(m_settings.tsProbesConfine.size() + "\tprobes have significant effects.");
//        reinit();
//        mapEQTLs();
//
//        // CONDITIONAL ANALYSES
//        // now determine top effect per probe..
//        // we'll regress out those effects later on.
//        runConditional(origOutputDir + "/SNP-Initial/", origOutputDir + "/Probe-Initial/", origOutputDir + "/SNP-Conditional/");
//        // now repeat the same procedure, but now conditioning for the gSNP effect
//        runConditional(origOutputDir + "/Probe-Initial/", origOutputDir + "/SNP-Initial/", origOutputDir + "/Probes-Conditional/");
//
//
//        // now summarize the hell out of the results...
//
//        // read initial snp probe pairs
//        eQTLObject[] snpInitial = readeQTLs(origOutputDir + "/SNP-Initial/eQTLsFDR.txt");
//        // read top effects for probes
//        eQTLObject[] probeInitial = readeQTLs(origOutputDir + "/Probe-Initial/eQTLsFDR.txt");
//        // read snp probe after conditioning
//        eQTLObject[] snpConditional = readeQTLs(origOutputDir + "/SNP-Conditional/eQTLsFDR.txt");
//        eQTLObject[] probeConditional = readeQTLs(origOutputDir + "/Probe-Conditional/eQTLsFDR.txt");
//
//        TextFile outfile = new TextFile(origOutputDir + "/ConditionalAnalysis.txt", TextFile.W);
//
//        String header = "probe\thugo\tgsnp\talleles\talleleassessed"
//                + "\tesnp\tesnpalleles\tesnpallelesasssesed"
//                + "\tgsnpp\tesnpp\tgsnpcondp\tesnpcondp"
//                + "\tgsnpmetab\tesnpmetab\tgsnpmetabcond\tesnpmetabcond"
//                + "\tgsnpb\tesnpb\tgsnpbcond\tesnpbcond"
//                + "\tgsnpss\tesnpss\tgsnpbconss\tesnpbconss"
//                + "\tgsnpsFDR\tesnpsFDR\tgsnpbconsFDR\tesnpbconsFDR";
//        outfile.writeln(header);
//
//        for (eQTLObject ei : snpInitial) {
//            // iterate through probeintitial
//            eQTLObject matchingpi = null;
//            eQTLObject matchingsc = null;
//            eQTLObject matchingpc = null;
//            for (eQTLObject pi : probeInitial) {
//                if (pi.probe.equals(ei.probe)) {
//                    // this is the thing to print
//                    matchingpi = pi;
//                    break;
//                }
//            }
//
//            // now the snp conditional on top probe
//            for (eQTLObject ec : snpConditional) {
//                if (ec.snp.equals(ei.snp) && ec.snp.equals(ei.snp)) {
//                    matchingsc = ec;
//                    break;
//                }
//            }
//
//            // now the probe conditional on top snp
//            for (eQTLObject pc : probeConditional) {
//                if (pc.snp.equals(matchingpi.snp) && pc.probe.equals(matchingpi.probe)) {
//                    matchingpc = pc;
//                    break;
//                }
//            }
//
//            String output = ei.probe
//                    + "\t" + ei.hugo
//                    + "\t" + ei.snp
//                    + "\t" + ei.alleles
//                    + "\t" + ei.alleleAssessed
//                    + "\t" + matchingpi.snp
//                    + "\t" + matchingpi.alleles
//                    + "\t" + matchingpi.alleleAssessed
//                    + "\t" + ei.p
//                    + "\t" + matchingpi.p
//                    + "\t" + matchingsc.p
//                    + "\t" + matchingpc.p
//                    + "\t" + ei.metab
//                    + "\t" + matchingpi.metab
//                    + "\t" + matchingsc.metab
//                    + "\t" + matchingpc.metab
//                    + "\t" + ei.b
//                    + "\t" + matchingpi.b
//                    + "\t" + matchingsc.b
//                    + "\t" + matchingpc.b
//                    + "\t" + ei.samplesize
//                    + "\t" + matchingpi.samplesize
//                    + "\t" + matchingsc.samplesize
//                    + "\t" + matchingpc.samplesize
//                    + "\t" + ei.FDR
//                    + "\t" + matchingpi.FDR
//                    + "\t" + matchingsc.FDR
//                    + "\t" + matchingpc.FDR;
//
//
//            outfile.writeln(output);
//
//
//        }
//        outfile.close();
	}
	
	private eQTLObject[] readeQTLs(String eqtlfile) throws IOException {
		ArrayList<eQTLObject> output = new ArrayList<eQTLObject>();
		TextFile tf = new TextFile(eqtlfile, TextFile.R);
		tf.readLine(); // skip header
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			
			eQTLObject e = new eQTLObject();
			e.p = elems[0];
			e.snp = elems[1];
			e.snpchr = elems[2];
			e.snpchrpos = elems[3];
			e.probe = elems[4];
			e.probechr = elems[5];
			e.probechrpos = elems[6];
			e.alleleAssessed = elems[QTLTextFile.ASESSEDALLELE];
			e.alleles = elems[QTLTextFile.ASESSEDALLELE - 1];
			e.samplesize = elems[QTLTextFile.DATASETSIZE];
			e.metab = elems[QTLTextFile.METAB];
			e.metaz = elems[QTLTextFile.METAZ];
			e.b = elems[QTLTextFile.DATASETB];
			e.hugo = elems[QTLTextFile.HUGO];
			e.FDR = elems[elems.length - 1];
			output.add(e);
			
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
//
//
		return output.toArray(new eQTLObject[0]);
	}
	
	private HashMap<Pair<String, String>, Double> readFDRFile(String file) throws IOException {
		HashMap<Pair<String, String>, Double> output = new HashMap<Pair<String, String>, Double>();
		TextFile tfSNPsCond = new TextFile(file, TextFile.R);
		tfSNPsCond.readLine();
		
		String[] elems = tfSNPsCond.readLineElems(TextFile.tab);
		while (elems != null) {
			String snp = elems[1];
			String probe = elems[4];
			String fdrStr = elems[elems.length - 1];
			Double fdr = Double.parseDouble(fdrStr);
			output.put(new Pair<String, String>(snp, probe), fdr);
			elems = tfSNPsCond.readLineElems(TextFile.tab);
		}
		tfSNPsCond.close();
		return output;
	}
	
	private HashMap<Pair<String, String>, EQTL> readEQTLFile(String eqtlfilename) throws IOException {
		HashMap<Pair<String, String>, EQTL> output = new HashMap<Pair<String, String>, EQTL>();
		TextFile efile = new TextFile(eqtlfilename, TextFile.R);
		efile.readLine();
		String[] elems = efile.readLineElems(TextFile.tab);
		while (elems != null) {
			
			EQTL e = new EQTL();
			
			e.setPvalue(Double.parseDouble(elems[0]));
			e.setRsName(elems[1]);
			e.setRsChr(ChrAnnotation.parseChr(elems[2]));
			e.setRsChrPos(Integer.parseInt(elems[3]));
			e.setProbe(elems[4]);
			e.setProbeChr(ChrAnnotation.parseChr(elems[5]));
			e.setProbeChrPos(Integer.parseInt(elems[6]));
			e.setAlleles(elems[8]);
			e.setAlleleAssessed(elems[9]);
			if (elems.length > 20) {
				if (!elems[20].equals("" + null)) {
					e.setFC(elems[20]);
				}
			}
			e.setZscore(Double.parseDouble(elems[QTLTextFile.METAZ]));
			e.setMetaBeta(elems[QTLTextFile.METAB]);
			e.setBeta(elems[QTLTextFile.DATASETB]);
			
			output.put(new Pair<String, String>(elems[1], elems[4]), e);
			elems = efile.readLineElems(TextFile.tab);
		}
		efile.close();
		
		return output;
	}
	
	void setSkipAllEQTLMapping() {
		skipAlleQTLMapping = true;
	}
	
	private class eQTLObject {
		
		String p;
		String snp;
		String probe;
		String metaz;
		String metab;
		String z;
		String b;
		String hugo;
		String samplesize;
		String alleleAssessed;
		String alleles;
		String probechr;
		String probechrpos;
		String snpchr;
		String snpchrpos;
		String FDR;
        /*
         Probe
         Probe Chr
         Probe ChrPos
         Probe HUGO
         gSNP
         gSNP Chr
         gSNP ChrPos
         gSNP alleles
         gSNP allele assessed
         eSNP
         eSNP Chr
         eSNP ChrPos
         eSNP alleles
         eSNP allele assessed
         LD
         gSNP-PVal
         eSNP-PVal
         gSNP-PVal-Conditional
         eSNP-PVal-Conditional
         gSNP-ZMeta
         eSNP-ZMeta
         gSNP-ZMeta-Conditional
         eSNP-ZMeta-Conditional
         gSNP-BMeta
         eSNP-BMeta
         gSNP-BMeta-Conditional
         eSNP-BMeta-Conditional
         gSNP-B
         eSNP-B
         gSNP-B-Conditional
         eSNP-B-Conditional
         gSNP-FC
         eSNP-FC
         gSNP-FC-Conditional
         eSNP-FC-Conditional
         gSNP-FDR
         eSNP-FDR
         gSNP-FDR-Conditional
         eSNP-FDR-Conditional 
         */
	}
	
	//    private void runConditional(String input, String condition, String output) throws IOException, Exception {
//        // now determine top effect per probe..
//        // we'll regress out those effects later on.
//        determineTopProbeEffect(condition);
//        // determine SNP-Probe combo's to test
//        HashMap<String, HashSet<String>> snpProbePairs = determineSNPProbeCombos(input + "/eQTLs.txt");
//        // run conditional analysis: test only SNP/Probe pairs, condition on top probe effect
//        m_settings.outputReportsDir = output;
//        m_settings.plotOutputDirectory = output + "/plots/";
//        Gpio.createDir(m_settings.plotOutputDirectory);
//        Gpio.createDir(m_settings.outputReportsDir);
//        m_settings.regressOutEQTLEffectFileName = condition + "/eQTLsTopProbeEffects.txt";
//        m_settings.tsSNPProbeCombinationsConfine = snpProbePairs;
//        m_settings.performEQTLAnalysisOnSNPProbeCombinationSubset = true;
//        reinit();
//        mapEQTLs();
//    }
	private void determineTopProbeEffect(String outdir) throws IOException {
		TextFile efile = new TextFile(outdir + "/eQTLs.txt", TextFile.R);
		TextFile efileout = new TextFile(outdir + "/eQTLsTopProbeEffects.txt", TextFile.W);
		efileout.writeln(efile.readLine());
		String[] elems = efile.readLineElems(TextFile.tab);
		HashSet<String> probesVisited = new HashSet<String>();
		while (elems != null) {
			if (!probesVisited.contains(elems[4])) {
				efileout.writeln(Strings.concat(elems, Strings.tab));
				probesVisited.add(elems[4]);
			}
			elems = efile.readLineElems(TextFile.tab);
		}
		efileout.close();
		efile.close();
	}
	
	private HashMap<String, HashSet<String>> determineSNPProbeCombos(String filename) throws IOException {
		HashMap<String, HashSet<String>> gsnpProbePairs = new HashMap<String, HashSet<String>>();
		TextFile efile = new TextFile(filename, TextFile.R);
		efile.readLine();
		String[] elems = efile.readLineElems(TextFile.tab);
		while (elems != null) {
			
			String snp = elems[1];
			String probe = elems[4];
			HashSet<String> probes = gsnpProbePairs.get(snp);
			if (probes == null) {
				probes = new HashSet<String>();
			}
			probes.add(probe);
			gsnpProbePairs.put(snp, probes);
			elems = efile.readLineElems(TextFile.tab);
		}
		efile.close();
		return gsnpProbePairs;
	}
}
