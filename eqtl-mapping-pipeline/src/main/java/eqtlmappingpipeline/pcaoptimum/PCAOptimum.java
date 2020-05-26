/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.pcaoptimum;

import eqtlmappingpipeline.metaqtl3.FDR;
import eqtlmappingpipeline.metaqtl3.MetaQTL3;
import eqtlmappingpipeline.metaqtl3.containers.Settings;
import eqtlmappingpipeline.normalization.Normalizer;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;
import umcg.genetica.console.ConsoleGUIElems;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDataset;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDatasetSettings;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.text.Strings;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * @author harmjan
 */
public class PCAOptimum extends MetaQTL3 {
	
	//	protected String inexpplatform;
//	protected String inexpannot;
//	protected String ingt;
//	protected String gte;
	protected Integer m_threads = 1;
	protected int permutations = 10;
	protected boolean covariatesremoved = false;
	protected String cissnps;
	protected String transsnps;
	private boolean performEigenVectorQTLMapping;
	
	public void setCovariatesRemoved(boolean b) {
		covariatesremoved = b;
	}
	
	public void setSNPSets(String cissnps, String transsnps) {
		this.cissnps = cissnps;
		this.transsnps = transsnps;
	}
	
	public void setPerformpcqtlNormalization(boolean performEigenvectorQTLMapping) {
		this.performEigenVectorQTLMapping = performEigenvectorQTLMapping;
	}
	
	@Override
	public void initialize(String xmlSettingsFile, String texttoreplace, String texttoreplacewith,
						   String ingt, String inexp, String inexpplatform, String inexpannot, String gte,
						   String out, boolean cis, boolean trans, int perm, boolean textout, boolean binout, String snpfile,
						   Integer threads, Integer maxNrResults, String regressouteqtls, String snpprobecombofile, boolean skipdotplot,
						   boolean skipqqplot, Long rseed, Double maf, Double hwe) throws IOException, Exception {
		
		HashSet<String> cisSnpsToTest = new HashSet<String>();
		HashSet<String> transSnpsToTest = new HashSet<String>();
		
		
		if (cissnps != null) {
			cis = true;
		} else {
			cis = false;
		}
		
		if (transsnps != null) {
			trans = true;
		} else {
			trans = false;
		}
		
		if (cis) {
			System.out.println("Loading cis SNP set from: " + cissnps);
			TextFile tf = new TextFile(cissnps, TextFile.R);
			cisSnpsToTest.addAll(tf.readAsArrayList());
			tf.close();
		}
		
		if (trans) {
			System.out.println("Loading trans SNP set from: " + transsnps);
			TextFile tf = new TextFile(transsnps, TextFile.R);
			transSnpsToTest.addAll(tf.readAsArrayList());
			tf.close();
		}
		
		if (xmlSettingsFile == null) {
			// store all settings in a settings object
			if (!out.endsWith("/")) {
				out += "/";
			}
			if (!Gpio.exists(out)) {
				Gpio.createDir(out);
			}
			
			permutations = perm;
			
			
			m_settings = new Settings();
			m_settings.numberOfVariantsToBuffer = 1000;
			if (snpfile != null) {
				m_settings.numberOfVariantsToBuffer = 1;
			}
			m_settings.fullFdrOutput = false;
			
			
			int nrProcs = Runtime.getRuntime().availableProcessors();
			if (threads != null && threads > 0 && threads <= nrProcs) {
				//
				m_threads = threads;
			} else {
				if (threads != null && threads > nrProcs) {
					System.out.println("The number of threads you set using the command line is not correct for your system. You set " + threads + " threads, while your machine has " + nrProcs + " processors");
				}
				threads = nrProcs;
			}
			
			m_threads = threads;
			int round = 0;
			
			//String nextInExp = origInExp + ".QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.txt.gz";
			String outputdir = out + "Cis-0PCAsRemoved";
			if (!outputdir.endsWith("/")) {
				outputdir += "/";
			}
			if (!Gpio.exists(outputdir)) {
				Gpio.createDir(outputdir);
			}
			
			m_settings.datasetSettings = new ArrayList<>();
			TriTyperGeneticalGenomicsDatasetSettings ds = new TriTyperGeneticalGenomicsDatasetSettings();
			ds.name = "Dataset";
			ds.probeannotation = inexpannot;
			ds.expressionplatform = inexpplatform;
			ds.expressionLocation = inexp;
			ds.genotypeLocation = ingt;
			ds.genotypeToExpressionCoupling = gte;
			m_settings.datasetSettings.add(ds);
			
		} else {
//			System.out.println("This tool does not support the use of XML settings files.");
//			System.exit(-1);
			m_settings = new Settings();
			m_settings.load(xmlSettingsFile);
			permutations = m_settings.nrPermutationsFDR;
			m_threads = m_settings.nrThreads;
			// store all settings in a settings object
			out = m_settings.outputReportsDir;
			if (!out.endsWith("/")) {
				out += "/";
			}
			if (!Gpio.exists(out)) {
				Gpio.createDir(out);
			}
		}
		
		ArrayList<String> origExpDs = new ArrayList<>();
		ArrayList<Integer> pcs = null;
		boolean alldshavesamepcs = true;
		
		for (int d = 0; d < m_settings.datasetSettings.size(); d++) {
			origExpDs.add(m_settings.datasetSettings.get(d).expressionLocation);
			if (d == 0) {
				pcs = getPCs(d);
			} else {
				ArrayList<Integer> dspcs = getPCs(d);
				HashSet<Integer> dspcshash = new HashSet<>();
				dspcshash.addAll(dspcs);
				for (Integer pc : pcs) {
					if (!(dspcs.contains(pc))) {
						System.out.println("Error: could not find file for " + pc + " for dataset " + m_settings.datasetSettings.get(d).name);
						alldshavesamepcs = false;
					}
				}
			}
		}
		
		
		if (!alldshavesamepcs) {
			System.out.println("Error: not all pc files present for all datasets. ");
			System.exit(-1);
		}
		
		Collections.sort(pcs);
		
		int max = 0;
		int stepSize = 0;
		for (int i = 0; i < (pcs.size() - 1); i++) {
			if (i == 0) {
				if (pcs.get(pcs.size() - 1) > max) {
					max = pcs.get(pcs.size() - 1);
				}
			}
			if (pcs.get(i) > max) {
				max = pcs.get(i);
			}
			stepSize += pcs.get(i + 1) - pcs.get(i);
		}
		
		if (pcs.size() > 1) {
			if ((((double) stepSize / (pcs.size() - 1)) % 1) != 0) {
				System.out.println("Step size is invalid."
						+ "\n Please look in to the input directory for missing files");
				System.out.println((((double) stepSize / (pcs.size() - 1)) % 1));
				System.exit(0);
			}
			stepSize = (int) ((double) stepSize / (pcs.size() - 1));
		} else {
			stepSize = pcs.get(0);
			max = pcs.get(0);
		}
		
		
		System.out.println("Determined max: " + max);
		System.out.println("Determined step size: " + stepSize);
		
		if (performEigenVectorQTLMapping) {
			performeQTLMappingOverEigenvectorMatrixAndReNormalize(out, stepSize, max, maxNrResults);
			
		}
		
		
		for (int pca = 0; pca <= max; pca += stepSize) {
			for (int d = 0; d < m_settings.datasetSettings.size(); d++) {
				String expfile = origExpDs.get(d);
				
				
				if (pca > 0) {
					String startExpressionFileName = expfile;
					File st = new File(startExpressionFileName);
					
					// strip the parent dir name
					String parentDir = Gpio.getParentDir(expfile);
					parentDir += Gpio.getFileSeparator();
					String minimalFilename = st.getName();
					if (minimalFilename.endsWith(".txt")) {
						minimalFilename = minimalFilename.substring(0, minimalFilename.length() - 4);
					} else if (minimalFilename.endsWith(".txt.gz")) {
						minimalFilename = minimalFilename.substring(0, minimalFilename.length() - 7);
					}
					if (performEigenVectorQTLMapping) {
						expfile = parentDir + minimalFilename + "." + pca + "PCAsOverSamplesRemoved-GeneticVectorsNotRemoved.txt.gz";
					} else {
						expfile = parentDir + minimalFilename + "." + pca + "PCAsOverSamplesRemoved.txt.gz";
					}
				}
				
				
				// check whether the file exists...
				if (!Gpio.exists(expfile)) {
					System.err.println("Could not find file for pca: " + pca + "\t" + expfile);
					System.exit(-1);
				}
				m_settings.datasetSettings.get(d).expressionLocation = expfile;
			}
			try {
				
				if (cis) {
					String outputDir = out + "Cis-" + pca + "PCAsRemoved/";
					if (performEigenVectorQTLMapping && pca > 0) {
						outputDir = out + "Cis-" + pca + "PCAsRemoved-GeneticVectorsNotRemoved/";
					}
					if ((pca == 0 && !Gpio.exists(outputDir + "eQTLProbesFDR0.05.txt.gz")) || pca > 0) {
						for (int d = 0; d < m_settings.datasetSettings.size(); d++) {
							System.out.println(m_settings.datasetSettings.get(d).toString());
							
						}
//						System.exit(-1);
						performeQTLMapping(true, false, outputDir, cisSnpsToTest, null, m_settings.nrThreads, m_settings.maxNrMostSignificantEQTLs);
						cleanup();
					}
				}
				if (trans) {
					String outputDir = out + "Trans-" + pca + "PCAsRemoved/";
					if (performEigenVectorQTLMapping && pca > 0) {
						outputDir = out + "Trans-" + pca + "PCAsRemoved-GeneticVectorsNotRemoved/";
					}
					if ((pca == 0 && !Gpio.exists(outputDir + "eQTLProbesFDR0.05.txt.gz")) || pca > 0) {
						performeQTLMapping(false, true, outputDir, transSnpsToTest, null, m_settings.nrThreads, m_settings.maxNrMostSignificantEQTLs);
						cleanup();
					}
				}
			} catch (NullPointerException e) {
				e.printStackTrace();
				System.exit(-1);
			}
		}
		
		
		PCAOptimumInventorize pi = new PCAOptimumInventorize();
		if (performEigenVectorQTLMapping) {
			pi.inventorypcqtl(out, cis, trans, max, stepSize);
		} else {
			pi.inventory(out, cis, trans, max, stepSize);
		}
		
	}
	
	private ArrayList<Integer> getPCs(int d) {
		String origInExp = m_settings.datasetSettings.get(d).expressionLocation;
		
		String parentDir = Gpio.getParentDir(origInExp);
		System.out.println("Looking for PCA corrected files in folder: " + parentDir);
		System.out.println("Base file name: " + origInExp);
		String[] fileList = Gpio.getListOfFiles(parentDir);
		ArrayList<Integer> dspcs = new ArrayList<Integer>();
		HashMap<Integer, String> stepToFile = new HashMap<Integer, String>();
		for (String f : fileList) {
			if (f.toLowerCase().contains("pcasoversamplesremoved") && !f.toLowerCase().contains("geneticvectorsnotremoved")) {
				String[] fileParts = f.split("\\.");
				for (String p : fileParts) {
					if (p.toLowerCase().contains("pcasoversamplesremoved")) {
						Integer pc = Integer.parseInt(p.toLowerCase().replace("pcasoversamplesremoved", ""));
						dspcs.add(pc);
						stepToFile.put(pc, f);
						System.out.println("Found file for PC: " + pc + "\t" + f);
						break;
					}
				}
			}
		}
		if (dspcs.isEmpty()) {
			System.out.println("No PCA corrected files."
					+ "\n Please first run the normalization procedure.");
			System.exit(0);
		}
		
		return dspcs;
	}
	
	protected void init() throws IOException, Exception {
		int numDatasets = m_settings.datasetSettings.size();
		m_gg = new TriTyperGeneticalGenomicsDataset[numDatasets];
		numAvailableInds = 0;
		for (int i = 0; i < numDatasets; i++) {
			
			System.out.println("- Loading dataset: " + m_settings.datasetSettings.get(i).name + "");
			System.out.println(ConsoleGUIElems.LINE);
			m_gg[i] = new TriTyperGeneticalGenomicsDataset(m_settings.datasetSettings.get(i));
			if (!m_settings.performParametricAnalysis) {
				m_gg[i].getExpressionData().rankAllExpressionData(m_settings.equalRankForTies);
			}
			m_gg[i].getExpressionData().calcAndSubtractMean();
			m_gg[i].getExpressionData().calcMeanAndVariance();
			
			numAvailableInds += m_gg[i].getExpressionToGenotypeIdArray().length;
			System.out.println(ConsoleGUIElems.LINE);
			System.out.println("");
		}
		
		
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
		
		
		// create threadss
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
		
		// save GTE's for future use
		for (int d = 0; d < m_gg.length; d++) {
			String outf = m_settings.outputReportsDir + "GTE-" + m_gg[d].getSettings().name + ".txt";
			TextFile tf = new TextFile(outf, TextFile.W);
			THashMap<String, String> samples = m_gg[d].getGenotypeToExpressionCouplings();
			for (Map.Entry<String, String> entry : samples.entrySet()) {
				tf.writeln(entry.getKey() + "\t" + entry.getValue());
			}
			tf.close();
		}
		
		printSummary();
	}
	
	protected void cleanup() {
		for (int i = 0; i < m_gg.length; i++) {
			this.m_gg[i] = null;
		}
		
		this.m_probeList = null;
		this.m_probeTranslationTable = null;
		this.m_snpList = null;
		this.m_snpTranslationTable = null;
		this.m_workPackages = null;
	}
	
	protected void performeQTLMapping(boolean cis, boolean trans, String out, HashSet<String> snpsToTest, THashSet<String> probesToTest, int threads, Integer maxNrResults) throws IOException, Exception {
//
		
		String outputdir = out;
		if (!Gpio.exists(outputdir)) {
			Gpio.createDir(outputdir);
		}
		// set output dir
		// set standard cis-settings
		Settings backup = m_settings;
		
		m_settings = new Settings();
		m_settings.datasetSettings = new ArrayList<>();
		for (int d = 0; d < backup.datasetSettings.size(); d++) {
			TriTyperGeneticalGenomicsDatasetSettings s = backup.datasetSettings.get(d);
			s.cisAnalysis = cis;
			s.transAnalysis = trans;
			if (probesToTest != null) {
				s.tsProbesConfine = probesToTest;
			}
			m_settings.datasetSettings.add(s);
			
		}
		
		
		m_settings.numberOfVariantsToBuffer = 1000;
		
		m_settings.createDotPlot = false;
		m_settings.displayWarnings = false;
		
		
		if (cis) {
			m_settings.ciseQTLAnalysMaxSNPProbeMidPointDistance = 250000;
		} else {
			m_settings.ciseQTLAnalysMaxSNPProbeMidPointDistance = 5000000;
		}
		
		m_settings.nrThreads = threads;
		m_settings.cisAnalysis = cis;
		m_settings.transAnalysis = trans;
		
		m_settings.nrPermutationsFDR = permutations;
		m_settings.tsSNPsConfine = snpsToTest;
		m_settings.numberOfVariantsToBuffer = 1000;
		if (snpsToTest != null) {
			m_settings.numberOfVariantsToBuffer = 1;
		}
		if (probesToTest != null) {
			m_settings.tsProbesConfine = probesToTest;
		}
		m_settings.outputReportsDir = outputdir;
		m_settings.createTEXTOutputFiles = true;
		m_settings.createBinaryOutputFiles = false;
		m_settings.randomNumberGenerator = new Random(m_settings.rSeed);
		m_settings.fdrType = FDR.FDRMethod.FULL;
		
		if (maxNrResults != null) {
			m_settings.maxNrMostSignificantEQTLs = maxNrResults;
			
		}
		
		init();
		// set standard trans settings
		super.mapEQTLs();
		cleanup();
		m_settings = backup;
	}
	
	protected void compareZScores(EQTL[] ciseqtls, EQTL[] transeqtls, EQTL[] originalCisEQTLs, EQTL[] originalTransEQTLs, String out, int pca) {
		HashMap<String, EQTL> origCis = new HashMap<String, EQTL>();
		HashMap<String, EQTL> origTrans = new HashMap<String, EQTL>();
		
		for (EQTL e : originalCisEQTLs) {
			origCis.put(e.getRsName() + "///" + e.getProbe(), e);
		}
		
		for (EQTL e : originalTransEQTLs) {
			origTrans.put(e.getRsName() + "///" + e.getProbe(), e);
		}
		
		int cisShared = 0;
		int transShared = 0;
		
		ArrayList<Double> cisXAL = new ArrayList<Double>();
		ArrayList<Double> cisYAL = new ArrayList<Double>();
		ArrayList<Double> transXAL = new ArrayList<Double>();
		ArrayList<Double> transYAL = new ArrayList<Double>();
		
		for (EQTL e : ciseqtls) {
			EQTL origE = origTrans.get(e.getRsName() + "///" + e.getProbe());
			if (origE != null) {
				cisShared++;
				cisXAL.add(e.getZscore());
				cisYAL.add(origE.getZscore());
			}
		}
		
		for (EQTL e : transeqtls) {
			EQTL origE = origTrans.get(e.getRsName() + "///" + e.getProbe());
			if (origE != null) {
				transShared++;
				transXAL.add(e.getZscore());
				transYAL.add(origE.getZscore());
			}
		}
		
		
		double[] cisx = toArray(cisXAL);
		double[] cisy = toArray(cisYAL);
		double[] transx = toArray(transXAL);
		double[] transy = toArray(transYAL);
		
		PCAOptimumPlot p = new PCAOptimumPlot(500, 1000, true, out);
		p.plot(cisx, cisy, transx, transy, originalCisEQTLs.length, ciseqtls.length, cisShared, originalTransEQTLs.length, transeqtls.length, transShared);
		p.draw(out);
		
	}
	
	protected double[] toArray(ArrayList<Double> x) {
		double[] y = new double[x.size()];
		int i = 0;
		for (Double z : x) {
			y[i] = z;
			i++;
		}
		return y;
	}
	
	public void performeQTLMappingOverEigenvectorMatrixAndReNormalize(String out, int stepSize, int max, Integer maxNrResults) throws IOException, Exception {
		Normalizer n = new Normalizer();
		
		ArrayList<String> origexp = new ArrayList<>();
		ArrayList<THashSet<String>> origprobeConfine = new ArrayList<>();
		for (int d = 0; d < m_settings.datasetSettings.size(); d++) {
			
			String origInExp = m_settings.datasetSettings.get(d).expressionLocation;
			origprobeConfine.add(m_settings.datasetSettings.get(d).tsProbesConfine);
			String parentDir = Gpio.getParentDir(origInExp);
			parentDir += Gpio.getFileSeparator();
			
			// check whether transposed eQTL file is there
			String expressionFileName = Gpio.getFileName(origInExp);
			if (parentDir == null) {
				parentDir = "";
			}
			
			if (expressionFileName.contains(".txt.gz")) {
				expressionFileName = expressionFileName.replaceAll(".txt.gz", "");
			} else {
				expressionFileName = expressionFileName.replaceAll(".txt", "");
			}
			
			origexp.add(origInExp);
			
			String outputFileNamePrefix = parentDir + expressionFileName;
			String eigenvectorFileTransposed = outputFileNamePrefix + ".PCAOverSamplesEigenvectorsTransposed.txt.gz";
			String eigenvectors = outputFileNamePrefix + ".PCAOverSamplesEigenvectors.txt.gz";
			System.out.println("Looking for: ");
			System.out.println(eigenvectors);
			System.out.println(eigenvectorFileTransposed);
			if (Gpio.exists(eigenvectorFileTransposed) && Gpio.exists(eigenvectorFileTransposed)) {
				System.out.println("Both files found in " + parentDir + ". No action required.");
			} else if (Gpio.exists(eigenvectors) && !Gpio.exists(eigenvectorFileTransposed)) {
				// transpose eigenvector matrix
				System.out.println("Eigenvector matrix found. Transposing " + eigenvectors);
				DoubleMatrixDataset<String, String> ds = DoubleMatrixDataset.loadDoubleData(eigenvectors);
				ds.saveDice(eigenvectorFileTransposed);
			} else {
				System.out.println("Files not found. Repeating PCA.");
				Triple<String, String, String> locations = n.calculatePcaOnly(origInExp);
				eigenvectorFileTransposed = locations.getLeft();
				String eigenvectorFile = locations.getMiddle();
				String pcaFile = locations.getRight();
			}
			m_settings.datasetSettings.get(d).expressionLocation = eigenvectorFileTransposed;
			
		}
		// Eigenvector mapping
		
		int nrToRemove = max + 1;
		
		THashSet<String> probesToTest = new THashSet<String>();
		for (int i = 1; i < nrToRemove; i++) {
			probesToTest.add("Comp" + i);
		}

//
		// ExpressionData.txt.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemoved.PCAOverSamplesEigenvectorsTransposed
		System.out.println("");
		System.out.println("Now performing QTL analysis on PCs");
		System.out.println("");
		performeQTLMapping(true, true, out + "CisTrans-PCAEigenVectors/", m_settings.tsSNPsConfine, probesToTest, m_threads, maxNrResults);
		cleanup();
		
		QTLTextFile etf = new QTLTextFile(out + "CisTrans-PCAEigenVectors/eQTLProbesFDR0.05.txt.gz", QTLTextFile.R);
		EQTL[] eigenvectorEQTLs = etf.read();
		etf.close();
		
		HashSet<Integer> geneticEigenVectors = new HashSet<Integer>();
		for (EQTL e : eigenvectorEQTLs) {
			Double fdr = e.getFDR();
			if (fdr != null) {
				if (fdr == 0) {
					String probe = e.getProbe();
					Integer compId = Integer.parseInt(probe.replace("Comp", ""));
					// quick hack: component 1 captures population stratification information...
					if (compId > 1) {
						geneticEigenVectors.add(compId);
					}
				}
			}
		}

//        System.out.println("Repeating PCA analysis, without removal of genetically controlled PCs");
		if (!geneticEigenVectors.isEmpty()) {
			System.out.println("These PCs are under genetic control: " + Strings.concat(geneticEigenVectors.toArray(new Integer[0]), Strings.comma));
			System.out.println();
		} else {
			System.out.println("No PCs are under genetic control.");
			System.out.println();
		}
		
		for (int d = 0; d < m_settings.datasetSettings.size(); d++) {
			String origInExp = origexp.get(d);
			String parentDir = Gpio.getParentDir(origInExp);
			parentDir = parentDir + Gpio.getFileSeparator();
			
			// check whether transposed eQTL file is there
			String expressionFileName = Gpio.getFileName(origInExp);
			if (parentDir == null) {
				parentDir = "";
			}
			
			if (expressionFileName.contains(".txt.gz")) {
				expressionFileName = expressionFileName.replaceAll(".txt.gz", "");
			} else {
				expressionFileName = expressionFileName.replaceAll(".txt", "");
			}
			
			
			String outputFileNamePrefix = parentDir + expressionFileName;
			String eigenvectorFileTranspose = outputFileNamePrefix + ".PCAOverSamplesEigenvectorsTransposed.txt.gz";
			String eigenvectorFile = outputFileNamePrefix + ".PCAOverSamplesEigenvectors.txt.gz";
			String pcaFile = outputFileNamePrefix + ".PCAOverSamplesPrincipalComponents.txt.gz";
			n.repeatPCAOmitCertainPCAs(geneticEigenVectors, parentDir, origInExp, eigenvectorFile, pcaFile, max, stepSize);
			m_settings.datasetSettings.get(d).expressionLocation = origInExp;
			m_settings.datasetSettings.get(d).tsProbesConfine = origprobeConfine.get(d);
		}
	}
	
	public void alternativeInitialize(String settingsfile, String ingt, String inexp, String inexpplatform, String inexpannot, String gte, String out,
									  boolean cis, boolean trans, int perm, String snpfile, Integer threads, boolean sortsnps) throws IOException, Exception {
		
		
		if (settingsfile != null) {
			m_settings = new Settings();
			m_settings.load(settingsfile);
			permutations = m_settings.nrPermutationsFDR;
			m_threads = m_settings.nrThreads;
			out = m_settings.outputReportsDir;
		} else {
			permutations = perm;
			m_settings = new Settings();
			m_settings.datasetSettings = new ArrayList<>();
			m_settings.sortsnps = sortsnps;
			m_settings.fullFdrOutput = false;
			m_settings.displayWarnings = false;
			int nrProcs = Runtime.getRuntime().availableProcessors();
			if (threads != null && threads > 0 && threads <= nrProcs) {
				//
				m_threads = threads;
			} else {
				if (threads != null && threads > nrProcs) {
					System.out.println("The number of threads you set using the command line is not correct for your system. You set " + threads + " threads, while your machine has " + nrProcs + " processors");
				}
				threads = nrProcs;
			}
			
			m_threads = threads;
			m_settings.createDotPlot = false;
			m_settings.createQQPlot = false;
			m_settings.fullFdrOutput = false;
		}
		
		if (!out.endsWith(Gpio.getFileSeparator())) {
			out += Gpio.getFileSeparator();
		}
		if (!Gpio.exists(out)) {
			Gpio.createDir(out);
		}
		if (cissnps != null) {
			cis = true;
		} else {
			cis = false;
		}
		
		if (transsnps != null) {
			trans = true;
		} else {
			trans = false;
		}
		
		m_settings.cisAnalysis = cis;
		m_settings.transAnalysis = trans;
		
		if (settingsfile == null) {
			TriTyperGeneticalGenomicsDatasetSettings ds = new TriTyperGeneticalGenomicsDatasetSettings();
			ds.expressionLocation = inexp;
			ds.probeannotation = inexpannot;
			ds.expressionplatform = inexpplatform;
			ds.genotypeLocation = ingt;
			ds.genotypeToExpressionCoupling = gte;
			ds.cisAnalysis = cis;
			ds.transAnalysis = trans;
		}
		
		if (snpfile != null) {
			TextFile f = new TextFile(snpfile, TextFile.R);
			m_settings.tsSNPsConfine = new HashSet<String>(f.readAsArrayList());
		}
		m_settings.numberOfVariantsToBuffer = 1000;
		if (snpfile != null) {
			m_settings.numberOfVariantsToBuffer = 1;
		}
		
	}
}
