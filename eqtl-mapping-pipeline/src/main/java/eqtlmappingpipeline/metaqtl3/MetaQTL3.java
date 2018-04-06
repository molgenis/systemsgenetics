/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl3;

import eqtlmappingpipeline.metaqtl3.containers.Settings;
import cern.colt.matrix.tint.IntMatrix2D;
import cern.colt.matrix.tint.impl.DenseIntMatrix2D;
import cern.colt.matrix.tint.impl.DenseLargeIntMatrix2D;
import com.itextpdf.text.DocumentException;
import eqtlmappingpipeline.metaqtl3.containers.WorkPackage;
import eqtlmappingpipeline.metaqtl3.containers.Result;
import umcg.genetica.math.stats.Descriptives;
import eqtlmappingpipeline.metaqtl3.graphics.EQTLDotPlot;
import eqtlmappingpipeline.metaqtl3.graphics.EQTLPlotter;

import java.io.IOException;
import java.util.*;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.logging.Level;
import java.util.logging.Logger;

import umcg.genetica.console.ConsoleGUIElems;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperExpressionData;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDataset;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDatasetSettings;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.util.RunTimer;

/**
 * @author harmjan
 */
public class MetaQTL3 {
	
	protected Settings m_settings;
	protected TriTyperGeneticalGenomicsDataset[] m_gg = null;
	protected String[] m_snpList;
	protected String[] m_probeList;
//    protected Integer[][] m_probeTranslationTable;
//    protected Integer[][] m_snpTranslationTable;
	
	//Defined -9 as null to not store it in an Integer
	protected IntMatrix2D m_probeTranslationTable;
	protected IntMatrix2D m_snpTranslationTable;
	
	protected int numAvailableInds;
	protected WorkPackage[] m_workPackages;
	private boolean dataHasCovariates;
	private Pair<List<String>, List<List<String>>> pathwayDefinitions;
	
	public MetaQTL3() {
	}
	
	public MetaQTL3(Settings settings) throws IOException, Exception {
		m_settings = settings;
		initialize(null, null, null, null, null, null, null, null, null, null, null, true, true, 0, true, false, null, null, null, null, null, false, false, null, null, null);
	}
	
	public void setOutputPlotThreshold(double d) {
		m_settings.plotOutputPValueCutOff = d;
		m_settings.plotOutputDirectory = m_settings.outputReportsDir;
		
	}
	
	public void initialize(String xmlSettingsFile, String texttoreplace, String texttoreplacewith, String texttoreplace2, String texttoreplace2with,
						   String ingt, String inexp, String inexpplatform, String inexpannot,
						   String gte, String out, boolean cis, boolean trans, int perm, boolean textout, boolean binout, String snpfile, Integer threads, Integer maxNrResults,
						   String regressouteqtls, String snpprobecombofile, boolean skipdotplot, boolean skipqqplot, Long rseed, Double maf, Double hwe) throws IOException, Exception {
		
		if (m_settings == null && xmlSettingsFile == null && ingt != null) {
			
			// check the settings
			boolean settingsOk = true;
			if (inexp == null || inexp.trim().length() == 0) {
				System.err.println("ERROR: you did not specify a gene expression file.");
				settingsOk = false;
			}
			
			if (inexpannot != null && inexpannot.trim().length() != 0) {
				if (inexpplatform == null || inexpplatform.trim().length() == 0) {
					System.err.println("ERROR: you specified " + inexpannot + " but you did not specify the platform (using --inexpplatform)!");
					settingsOk = false;
				}
			}
			
			if (out == null || out.trim().length() == 0) {
				System.err.println("ERROR: you did not specify an output directory.");
				settingsOk = false;
			}
			
			if (!settingsOk) {
				System.out.println();
				System.exit(0);
			}
			
			m_settings = new Settings();
			
			TriTyperGeneticalGenomicsDatasetSettings s = new TriTyperGeneticalGenomicsDatasetSettings();
			
			s.name = "Dataset";
			s.expressionLocation = inexp;
			s.expressionplatform = inexpplatform;
			s.probeannotation = inexpannot;
			s.genotypeLocation = ingt;
			s.genotypeToExpressionCoupling = gte;
			s.cisAnalysis = cis;
			s.transAnalysis = trans;
			
			m_settings.cisAnalysis = cis;
			m_settings.transAnalysis = trans;
			
			boolean cistrans = false;
			if (m_settings.cisAnalysis && m_settings.transAnalysis) {
				m_settings.confineProbesToProbesMappingToAnyChromosome = true;
			}
			
			m_settings.datasetSettings = new ArrayList<TriTyperGeneticalGenomicsDatasetSettings>();
			
			m_settings.regressOutEQTLEffectFileName = regressouteqtls;
			m_settings.datasetSettings.add(s);
			m_settings.nrThreads = threads;
			m_settings.cisAnalysis = cis;
			m_settings.transAnalysis = trans;
			m_settings.nrPermutationsFDR = perm;
			if (maf != null && maf > 0.0d && maf < 1.0d) {
				m_settings.snpQCMAFThreshold = maf;
			}
			if (hwe != null && hwe >= 0.0d && hwe <= 1.0d) {
				m_settings.snpQCHWEThreshold = hwe;
			}
			if (!out.endsWith("/")) {
				out += "/";
			}
			if (!Gpio.exists(out)) {
				Gpio.createDir(out);
			}
			if (snpfile != null) {
				m_settings.tsSNPsConfine = new HashSet<String>();
				TextFile ts = new TextFile(snpfile, TextFile.R);
				m_settings.strConfineSNP = snpfile;
				m_settings.tsSNPsConfine.addAll(ts.readAsArrayList());
				ts.close();
			}
			
			if (snpprobecombofile != null) {
				m_settings.loadSNPProbeConfinement(snpprobecombofile);
			}
			
			m_settings.outputReportsDir = out;
			m_settings.createTEXTOutputFiles = textout;
			m_settings.createBinaryOutputFiles = binout;
			if (maxNrResults != null && maxNrResults > 0) {
				m_settings.maxNrMostSignificantEQTLs = maxNrResults;
			}
			
			m_settings.createDotPlot = !skipdotplot;
			m_settings.createQQPlot = !skipqqplot;
			
			if (rseed != null) {
				m_settings.rSeed = rseed;
				m_settings.randomNumberGenerator = new Random(m_settings.rSeed);
			}
			
		} else if (m_settings == null && xmlSettingsFile != null) {
			// parse settings
			m_settings = new Settings();
			m_settings.settingsTextReplaceWith = texttoreplacewith;
			m_settings.settingsTextToReplace = texttoreplace;
			m_settings.settingsTextReplace2With = texttoreplace2with;
			m_settings.settingsTextToReplace2 = texttoreplace2;
			m_settings.load(xmlSettingsFile);
		} else if (m_settings == null) {
			System.out.println("ERROR: No input specified");
			System.exit(0);
		}
		
		// initialize dataset
		if (!m_settings.cisAnalysis && !m_settings.transAnalysis) {
			System.err.println("! WARNING: defaulting to CIS analysis (override with --trans or --trans and --cis))");
			m_settings.cisAnalysis = true;
		}

//        m_settings.randomNumberGenerator = new Random(m_settings.rSeed);
		m_settings.writeSettingsToDisk();
		
		int numDatasets = m_settings.datasetSettings.size();
		m_gg = new TriTyperGeneticalGenomicsDataset[numDatasets];
		numAvailableInds = 0;
		int nrOfDatasetsWithGeneExpressionData = 0;
		int nrDatasetsWithCovariates = 0;
		
		for (int i = 0; i < numDatasets; i++) {
			String covariateFile = m_settings.datasetSettings.get(i).covariateFile;
			if (covariateFile != null && Gpio.exists(covariateFile)) {
				nrDatasetsWithCovariates++;
			}
		}
		if (nrDatasetsWithCovariates >= 1 && nrDatasetsWithCovariates != m_gg.length) {
			System.err.println("Covariate files have not been specified for all datasets.");
			System.exit(-1);
		}
		
		if (nrDatasetsWithCovariates >= 1) {
			dataHasCovariates = true;
		}
		
		List<String> pathwayNames;
		List<List<String>> ensgsInPathways;
		if (m_settings.pathwayDefinition != null) {
			pathwayNames = new ArrayList<String>();
			ensgsInPathways = new ArrayList<List<String>>();
			if (Gpio.exists(m_settings.pathwayDefinition)) {
				TextFile tf = new TextFile(m_settings.pathwayDefinition, TextFile.R);
				String line;
				while ((line = tf.readLine()) != null) {
					List<String> ensgsThisPathway = new ArrayList<String>();
					String[] split = line.split("\t");
					for (int i = 2; i < split.length; i++) {
						ensgsThisPathway.add(split[i]);
					}
					pathwayNames.add(split[0]);
					ensgsInPathways.add(ensgsThisPathway);
				}
				tf.close();
				System.out.println("Read " + pathwayNames.size() + " pathways from " + m_settings.pathwayDefinition);
				pathwayDefinitions = new Pair<List<String>, List<List<String>>>(pathwayNames, ensgsInPathways);
				
				// be sure this is a cis-trans analysis on pathways...
				m_settings.cisAnalysis = true;
				m_settings.transAnalysis = true;
				for (int d = 0; d < numDatasets; d++) {
					// hooray for redundant settings: need to fix at some point.
					m_settings.datasetSettings.get(d).cisAnalysis = true;
					m_settings.datasetSettings.get(d).transAnalysis = true;
				}
			} else {
				System.err.println("Pathway defnition defined as: " + m_settings.pathwayDefinition + ", but file does not exist.");
				System.exit(-1);
			}
		} else {
			pathwayDefinitions = null;
		}
		
		
		for (int i = 0; i < numDatasets; i++) {
			System.out.println("- Loading dataset: " + m_settings.datasetSettings.get(i).name + "");
			m_settings.datasetSettings.get(i).confineProbesToProbesMappingToAnyChromosome = m_settings.confineProbesToProbesMappingToAnyChromosome;
			System.out.println(ConsoleGUIElems.LINE);
			m_gg[i] = new TriTyperGeneticalGenomicsDataset(m_settings.datasetSettings.get(i), pathwayDefinitions, m_settings.displayWarnings);
			
			if (m_gg[i].isExpressionDataLoadedCorrectly()) {
				nrOfDatasetsWithGeneExpressionData++;
			}
		}
		
		if (nrOfDatasetsWithGeneExpressionData == 0) {
			System.out.println("Error: none of your datasets contain any gene expression data for the settings you have specified");
			System.exit(0);
		}
		
		if (nrOfDatasetsWithGeneExpressionData != m_gg.length) {
			System.out.println("WARNING: was able to load gene expression data for " + nrOfDatasetsWithGeneExpressionData + " while you specified " + m_gg.length + " datasets in the settings.");
			
			// remove the datasets without expression data.
			TriTyperGeneticalGenomicsDataset[] tmp_gg = new TriTyperGeneticalGenomicsDataset[nrOfDatasetsWithGeneExpressionData];
			int ctr = 0;
			for (TriTyperGeneticalGenomicsDataset d : m_gg) {
				if (d.isExpressionDataLoadedCorrectly()) {
					tmp_gg[ctr] = d;
				}
			}
			m_gg = tmp_gg;
		}
		
		// rank and normalize data
		for (int i = 0; i < numDatasets; i++) {
			if (!m_settings.performParametricAnalysis) {
				m_gg[i].getExpressionData().rankAllExpressionData(m_settings.equalRankForTies);
			}
			m_gg[i].getExpressionData().calcAndSubtractMean();
			m_gg[i].getExpressionData().calcMeanAndVariance();
			numAvailableInds += m_gg[i].getExpressionToGenotypeIdArray().length;
		}
		
		if (m_settings.regressOutEQTLEffectFileName != null && m_settings.regressOutEQTLEffectFileName.trim().length() > 0) {
			if (!Gpio.exists(m_settings.regressOutEQTLEffectFileName)) {
				System.err.println("ERROR: you have specified an eQTL file to regress out, but the file was not found " + m_settings.regressOutEQTLEffectFileName);
				System.exit(0);
			}
			EQTLRegression eqr = new EQTLRegression();
			eqr.regressOutEQTLEffects(m_settings.regressOutEQTLEffectFileName, m_settings.regressOutEQTLEffectsSaveOutput, m_gg);
			
			// shouldn't we re-rank?
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
		
		
		System.out.println("Accumulating available data...");
		System.out.print(ConsoleGUIElems.LINE);
		
		createSNPList();
		createProbeList();
		
		// create WorkPackage objects
		determineSNPProbeCombinations();
//
//        if (m_settings.maxNrMostSignificantEQTLs != maxNrOfTests) {
//            m_settings.maxNrMostSignificantEQTLs = (int) maxNrOfTests;
//        }
		
		if (m_workPackages == null || m_workPackages.length == 0) {
			System.err.println("Error: No work detected");
			System.exit(0);
		}

//		// determine number of threads
		if (m_settings.nrThreads == null) {
			m_settings.nrThreads = Runtime.getRuntime().availableProcessors();
		}
//		else {
//			int numProcs = Runtime.getRuntime().availableProcessors();
//			if (m_settings.nrThreads > numProcs || m_settings.nrThreads < 1) {
//				m_settings.nrThreads = numProcs;
//			}
//		}
		
		if (m_workPackages.length < m_settings.nrThreads) {
			m_settings.nrThreads = m_workPackages.length;
		}
		printSummary();
		
	}
	
	protected void createSNPList() throws IOException {
		ArrayList<String> availableSNPs = new ArrayList<String>();
		HashSet<String> chrYSNPs = new HashSet<String>();
		HashSet<String> unknownPos = new HashSet<String>();
		HashSet<String> unknownchr = new HashSet<String>();
		
		TextFile excludedSNPs = new TextFile(m_settings.outputReportsDir + "excludedSNPsBySNPFilter.txt.gz", TextFile.W);
		
		HashSet<String> tmpAvailableSNPs = new HashSet<String>();
		HashSet<String> duplicatesWithinDataset = new HashSet<String>();
		for (TriTyperGeneticalGenomicsDataset m_gg1 : m_gg) {
			String[] snps = m_gg1.getGenotypeData().getSNPs();
			HashSet<String> seenSNPs = new HashSet<String>();
			for (String s : snps) {
				if (m_settings.tsSNPsConfine == null || m_settings.tsSNPsConfine.contains(s)) {
					if (seenSNPs.contains(s)) {
						duplicatesWithinDataset.add(s);
					} else {
						seenSNPs.add(s);
					}
					tmpAvailableSNPs.add(s);
				}
			}
		}
		
		String[] snps = tmpAvailableSNPs.toArray(new String[0]);
		for (int s = 0; s < snps.length; s++) {
			String snpName = snps[s];
			
			if (m_settings.tsSNPsConfine == null || m_settings.tsSNPsConfine.contains(snpName)) {
				StringBuilder reason = new StringBuilder(snps[s]).append("\t");
				if (m_gg.length == 1) {
					TriTyperGenotypeData ds = m_gg[0].getGenotypeData();
					Integer snpId = ds.getSnpToSNPId().get(snpName);
					Byte chr1 = ds.getChr(snpId);
					int chrPos1 = ds.getChrPos(snpId);
					boolean excludeSNP = false;
					
					if (chr1 >= 24) {
						chrYSNPs.add(snpName);
						reason.append("\tSNP is located on Y, MT, XY chromosome");
						excludeSNP = true;
					} else if (chrPos1 < 0) {
						unknownPos.add(snpName);
						reason.append("\tSNP has unknown mapping");
						excludeSNP = true;
					}
					
					if (duplicatesWithinDataset.contains(snpName)) {
						reason.append("\tSNP is present twice in one of the datasets");
						excludeSNP = true;
					}
					
					if (excludeSNP) {
						excludedSNPs.writeln(reason.toString());
					} else {
						availableSNPs.add(snpName);
					}
				} else { // meta-analysis requires a bit more extensive comparison
					boolean identicalMapping = true;
					boolean presentInAllDatasets = true;
					Byte chr1 = null;
					int chrPos1 = -1;
					
					boolean excludeSNP = false;

// check for identical mapping
					for (int d = 0; d < m_gg.length; d++) {
						TriTyperGenotypeData ds = m_gg[d].getGenotypeData();
						Integer snpId = ds.getSnpToSNPId().get(snpName);
						if (snpId == -9) {
							presentInAllDatasets = false;
							reason.append(";Not present in dataset ").append(d);
						} else {
							Byte chr2 = ds.getChr(snpId);
							int chrPos2 = ds.getChrPos(snpId);
							if (chr1 == null && chr2 != null) {
								chr1 = ds.getChr(snpId);
								chrPos1 = ds.getChrPos(snpId);
							} else {
								if (chr1 != null && chr2 == null) {
									identicalMapping = false;
									reason.append(";SNP has no chromosome position in dataset: ").append(d);
								} else if (chr1 != null && chr2 != null) {
									// compare positions
									if (!chr1.equals(chr2)) {
										identicalMapping = false;
										reason.append(";SNP maps to different chromosome in dataset: ").append(d).append("-chr").append(chr2);
									} else {
										if (chrPos1 != chrPos2) {
											identicalMapping = false;
											reason.append(";SNP maps to different chromosome position in dataset: ").append(d).append("-chr").append(chr2);
										}
									}
								}
							}
							
						}
					}
					
					// process the snps, based on the settings
					if (m_settings.confineToSNPsThatMapToChromosome != null) {
						if (chr1 == null || !chr1.equals(m_settings.confineToSNPsThatMapToChromosome)) {
							reason.append("\tSNP does not map to chromosome ").append(m_settings.confineToSNPsThatMapToChromosome).append(": chr").append(chr1);
							excludeSNP = true;
						}
					}
					if (m_settings.confineSNPsToSNPsPresentInAllDatasets && !presentInAllDatasets) {
						reason.append("\tSNP is not present in all datasets.");
						excludeSNP = true;
					}
					if (chr1 == null) {
						unknownchr.add(snpName);
						reason.append("\tSNP has unknown chromosome");
						excludeSNP = true;
					} else if (chr1 >= 24) {
						chrYSNPs.add(snpName);
						reason.append("\tSNP is located on Y, MT, XY chromosome");
						excludeSNP = true;
					} else if (chrPos1 < 0) {
						unknownPos.add(snpName);
						reason.append("\tSNP has unknown mapping");
						excludeSNP = true;
					}
					if (!identicalMapping) {
						reason.append("\tSNP has different mapping in different datasets");
						excludeSNP = true;
					}
					
					if (duplicatesWithinDataset.contains(snpName)) {
						reason.append("\tSNP is present twice in one of the datasets");
						excludeSNP = true;
					}
					
					if (excludeSNP) {
						excludedSNPs.writeln(reason.toString());
					} else {
						availableSNPs.add(snpName);
					}
				}
				
			}
			
		}
		
		System.out.println("- " + chrYSNPs.size() + " chromosome Y, MT or XY SNPs, " + unknownPos.size() + " SNPS with unknown position, " + unknownchr.size() + " with unknown chromosome.");
		System.out.println("- Remaining SNPs: " + availableSNPs.size());
		
		m_snpList = availableSNPs.toArray(new String[0]);
		
		
		// can we sort these based on position?
		class InSNP implements Comparable<InSNP> {
			
			String snp;
			int chr;
			int pos;
			
			public InSNP(String s, int chr, int pos) {
				this.snp = s;
				this.chr = chr;
				this.pos = pos;
			}
			
			@Override
			public int compareTo(InSNP o) {
				if (this.equals(o)) {
					return 0;
				} else if (this.chr > o.chr) {
					return 1;
				} else if (this.chr < o.chr) {
					return -1;
				} else {
					if (this.pos > o.pos) {
						return 1;
					} else if (this.pos < o.pos) {
						return -1;
					}
				}
				return 0;
			}
			
			@Override
			public boolean equals(Object o) {
				if (this == o) return true;
				if (o == null || getClass() != o.getClass()) return false;
				
				InSNP inSNP = (InSNP) o;
				
				if (chr != inSNP.chr) return false;
				if (pos != inSNP.pos) return false;
				return snp != null ? snp.equals(inSNP.snp) : inSNP.snp == null;
			}
			
			@Override
			public int hashCode() {
				int result = snp != null ? snp.hashCode() : 0;
				result = 31 * result + chr;
				result = 31 * result + pos;
				return result;
			}
		}
		
		System.out.println(m_snpList.length + " snps before sorting");
		ArrayList<InSNP> snpsArr = new ArrayList<>();
		for (int s = 0; s < m_snpList.length; s++) {
			int pos = 0;
			int chr = 0;
			String snp = m_snpList[s];
			for (int d = 0; d < m_gg.length; d++) {
				Integer id = m_gg[d].getGenotypeData().getSnpToSNPId().get(snp);
				if (id >= 0) {
					chr = m_gg[d].getGenotypeData().getChr(id);
					pos = m_gg[d].getGenotypeData().getChrPos(id);
				}
			}
			snpsArr.add(new InSNP(snp, chr, pos));
		}
		Collections.sort(snpsArr);
		System.out.println(snpsArr.size() + " snps after sorting");
		m_snpList = new String[m_snpList.length];
		for (int i = 0; i < m_snpList.length; i++) {
			m_snpList[i] = snpsArr.get(i).snp;
		}
		
		
		// create snp translation table..
		if ((m_gg.length * (long) m_snpList.length) < (Integer.MAX_VALUE - 2)) {
			m_snpTranslationTable = new DenseIntMatrix2D(m_gg.length, m_snpList.length);
		} else {
			m_snpTranslationTable = new DenseLargeIntMatrix2D(m_gg.length, m_snpList.length);
		}

//        m_snpTranslationTable = new Integer[m_gg.length][m_snpList.length];
		for (int p = 0; p < m_snpList.length; p++) {
			String snp = m_snpList[p];
			for (int d = 0; d < m_gg.length; d++) {
				Integer tmp = m_gg[d].getGenotypeData().getSnpToSNPId().get(snp);
				if (tmp == -9) {
					m_snpTranslationTable.setQuick(d, p, -9);
				} else {
					m_snpTranslationTable.setQuick(d, p, tmp);
				}
				
			}
		}
		excludedSNPs.close();
		
		// now determine which of the SNPs that was queried for does not exist in any of the datasets.
		if (m_settings.tsSNPsConfine != null) {
			Iterator<String> it = m_settings.tsSNPsConfine.iterator();
			if (m_settings.tsSNPsConfine.isEmpty()) {
				System.err.println("ERROR: a SNP confinement file is specified in the settings, but it is apparently empty? " + m_settings.strConfineSNP);
			} else {
				String next = it.next();
				
				TextFile querySNPNotPresent = new TextFile(m_settings.outputReportsDir + "querySNPsNotPresentInDataset.txt.gz", TextFile.W);
				while (it.hasNext()) {
					boolean isPresentInAnyDataset = false;
					for (int d = 0; d < m_gg.length; d++) {
						Integer id = m_gg[d].getGenotypeData().getSnpToSNPId().get(next);
						if (id != -9) {
							isPresentInAnyDataset = true;
						}
					}
					if (!isPresentInAnyDataset) {
						querySNPNotPresent.writeln(next);
					}
					next = it.next();
				}
				querySNPNotPresent.close();
			}
		}
		
	}
	
	protected void createProbeList() throws IOException {
		TextFile probeLog = new TextFile(m_settings.outputReportsDir + "ProbeQCLog.txt.gz", TextFile.W);
		
		System.out.println("- Determining available probes.");
		
		System.out.println("\t- Saving logfile to: " + m_settings.outputReportsDir + "ProbeQCLog.txt.gz");
		if (m_settings.confineProbesToProbesPresentInAllDatasets) {
			System.out.println("\t- Confining to probes present in all datasets");
		} else {
			System.out.println("\t- Not confining to probes present in all datasets");
		}
		
		if (m_settings.confineProbesToProbesMappingToAnyChromosome) {
			System.out.println("\t- Confining to probes that map to any chromosome (including probes without a valid position)");
		} else {
			System.out.println("\t- Confining to probes that map to autosomes, X, Y and MT chromosomes");
		}
		
		if (m_settings.confineToProbesThatMapToChromosome != null) {
			System.out.println("\t- Confining to probes that map to chromosome: " + m_settings.confineToProbesThatMapToChromosome);
		}
		
		int numExcluded = 0;
		// first check all datasets for duplicate probes...
		
		HashSet<String> duplicateProbes = new HashSet<String>();
		HashSet<String> allAvailableProbes = new HashSet<String>();
		for (int d = 0; d < m_gg.length; d++) {
			HashSet<String> availableProbes = new HashSet<String>();
			String[] probes = m_gg[d].getExpressionData().getProbes();
			for (String probe : probes) {
				if (availableProbes.contains(probe)) {
					duplicateProbes.add(probe);
					probeLog.writeln("Removing probe:\t" + probe + "\tis a duplicate");
				}
				availableProbes.add(probe);
				allAvailableProbes.add(probe);
			}
		}
		
		System.out.println("\t- " + allAvailableProbes.size() + " available probes for all datasets. Will now look for duplicate probes.");
		// remove duplicates from the analysis...
		
		String[] duplicates = duplicateProbes.toArray(new String[duplicateProbes.size()]);
		for (String duplicate : duplicates) {
			allAvailableProbes.remove(duplicate);
		}
		System.out.println("\t- " + allAvailableProbes.size() + " available probes for all datasets after removing " + duplicateProbes.size() + " duplicates");
		
		// determine in how many datasets the probes are present
		String[] availableProbeArray = allAvailableProbes.toArray(new String[0]);
		
		if (m_settings.tsProbesConfine == null) {
		} else {
			System.out.println("Probe confinement list has " + m_settings.tsProbesConfine.size() + " probes");
			availableProbeArray = m_settings.tsProbesConfine.toArray(new String[0]);
		}
		
		int mappingToDifferentPositionsAcrossDatasets = 0;
		int mapToWrongChromosome = 0;
		int invalidMappingPosition = 0;
		int nrProbesNotInAllDatasets = 0;
		HashSet<String> finalProbeList = new HashSet<String>();
		for (String probe : availableProbeArray) {
			int presence = 0;
			
			Byte chr = null;
			int chrPosStart = -1;
			int chrPosEnd = -1;
			
			boolean hasIdenticalMappingAcrossDatasets = true;
			String mappingOutput = "";
			for (int d = 0; d < m_gg.length; d++) {
				Integer probeId = m_gg[d].getExpressionData().getProbeToId().get(probe);
				if (probeId != -9) {
					presence++;
					if (chr == null) {
						chr = m_gg[d].getExpressionData().getChr()[probeId];
						if (chr == -1) {
							chr = null;
						}
						chrPosStart = m_gg[d].getExpressionData().getChrStart()[probeId];
						chrPosEnd = m_gg[d].getExpressionData().getChrStop()[probeId];
						
						mappingOutput += "\t" + m_gg[d].getSettings().name + ": Chr " + chr + "; Pos " + chrPosStart + "-" + chrPosEnd;
						
					} else {
						
						Byte chr2 = m_gg[d].getExpressionData().getChr()[probeId];
						if (chr2 == -1) {
							chr2 = null;
						}
						int chrPosStart2 = m_gg[d].getExpressionData().getChrStart()[probeId];
						int chrPosEnd2 = m_gg[d].getExpressionData().getChrStop()[probeId];
						
						if (chr2 == null) {
							hasIdenticalMappingAcrossDatasets = false;
						} else if (chrPosStart2 == -1 && chrPosStart != -1) {
							hasIdenticalMappingAcrossDatasets = false;
						} else if (chrPosEnd2 == -1 && chrPosEnd != -1) {
							hasIdenticalMappingAcrossDatasets = false;
						} else if (!chr.equals(chr2)) {
							hasIdenticalMappingAcrossDatasets = false;
						} else if (chrPosStart != (chrPosStart2)) {
							hasIdenticalMappingAcrossDatasets = false;
						} else if (chrPosEnd != (chrPosEnd2)) {
							hasIdenticalMappingAcrossDatasets = false;
						}
						mappingOutput += "\t" + m_gg[d].getSettings().name + ": Chr " + chr2 + "; Pos " + chrPosStart2 + "-" + chrPosEnd2;
					}
				}
			}
			
			// check whether the probe has a proper mapping
			boolean includeProbe = true;
			// check whether we would want to include this probe on the basis of its presence in other datasets
//	    if (m_settings.cisAnalysis && m_settings.transAnalysis && m_settings.confineProbesToProbesPresentInAllDatasets) { // ?
			
			if (m_settings.tsProbesConfine != null && m_settings.tsProbesConfine.contains(probe)) {
				includeProbe = true;
			} else if (m_settings.cisAnalysis && m_settings.transAnalysis) { // ?
				includeProbe = true;
			} else if (m_settings.confineProbesToProbesPresentInAllDatasets && presence < m_gg.length) {
				includeProbe = false;
				nrProbesNotInAllDatasets++;
				probeLog.writeln("Removing probe:\t" + probe + "\tpresent in " + presence + " / " + m_gg.length + "\tdatasets");
			} else if (chr == null || chr >= 25 || chr <= 0) {
				if (!m_settings.confineProbesToProbesMappingToAnyChromosome) {
					includeProbe = false;
					invalidMappingPosition++;
					probeLog.writeln("Removing probe:\t" + probe + "\t has no valid mapping position in any dataset: " + mappingOutput);
				}
			} else if (m_settings.confineToProbesThatMapToChromosome != null && !chr.equals(m_settings.confineToProbesThatMapToChromosome)) {
				// check whether this chromosome was requested to be analysed
				includeProbe = false;
				mapToWrongChromosome++;
				probeLog.writeln("Removing probe:\t" + probe + "\tmaps to wrong chromosome: " + mappingOutput);
			} else if (!hasIdenticalMappingAcrossDatasets) {
				// exclude the probe if it does not have identical mappings across datasets
				includeProbe = false;
				mappingToDifferentPositionsAcrossDatasets++;
				probeLog.writeln("Removing probe:\t" + probe + "\tmaps to different positions in datasets: " + mappingOutput);
			}
			
			if (includeProbe) {
				finalProbeList.add(probe);
			}
		}
		
		System.out.println("\t- " + finalProbeList.size() + "\tprobes finally included: ");
		if (m_settings.confineProbesToProbesPresentInAllDatasets) {
			System.out.println("\t\t- " + nrProbesNotInAllDatasets + " are not present in all datasets");
		}
		if (!m_settings.confineProbesToProbesMappingToAnyChromosome) {
			System.out.println("\t\t- " + invalidMappingPosition + " have an invalid mapping position");
		}
		System.out.println("\t\t- " + mappingToDifferentPositionsAcrossDatasets + " probes map to different positions in one or more datasets");
		if (m_settings.confineToProbesThatMapToChromosome != null) {
			System.out.println("\t\t- " + mapToWrongChromosome + " probes map to a different chromosome than the one selected (Chr: " + m_settings.confineToProbesThatMapToChromosome + ")");
		}
		
		if (finalProbeList.isEmpty()) {
			System.err.println("Error: no probes remaining after filter. Are your settings correct?");
			probeLog.close();
			System.exit(0);
		}
		m_probeList = finalProbeList.toArray(new String[finalProbeList.size()]);
		
		// create probe translation table..
		if ((m_gg.length * (long) m_probeList.length) < (Integer.MAX_VALUE - 2)) {
			m_probeTranslationTable = new DenseIntMatrix2D(m_gg.length, m_probeList.length);
		} else {
			m_probeTranslationTable = new DenseLargeIntMatrix2D(m_gg.length, m_probeList.length);
		}
		
		for (int p = 0; p < m_probeList.length; p++) {
			String probe = m_probeList[p];
			for (int d = 0; d < m_gg.length; d++) {
				Integer tmp = m_gg[d].getExpressionData().getProbeToId().get(probe);
				if (tmp == -9) {
					m_probeTranslationTable.setQuick(d, p, -9);
				} else {
					m_probeTranslationTable.setQuick(d, p, tmp);
				}
			}
		}
		
		probeLog.close();
	}
	
	public void mapEQTLs() throws IOException {
		
		// create work packages
		RunTimer t = new RunTimer();
		
		CalculationThread[] pool = new CalculationThread[m_settings.nrThreads];
		
		SNPLoader[] snploaders = new SNPLoader[m_gg.length];
		TriTyperExpressionData[] expressiondata = new TriTyperExpressionData[m_gg.length];
		for (int d = 0; d < snploaders.length; d++) {
			snploaders[d] = m_gg[d].getGenotypeData().createSNPLoader(m_settings.numberOfVariantsToBuffer);
			expressiondata[d] = m_gg[d].getExpressionData();
		}
		
		// initialize lookup tables
		int maxNrSamples = 0;
		for (int d = 0; d < m_gg.length; d++) {
			if (m_gg[d].getExpressionToGenotypeIdArray().length > maxNrSamples) {
				maxNrSamples = m_gg[d].getExpressionToGenotypeIdArray().length;
			}
		}
		Correlation.correlationToZScore(maxNrSamples);
		Descriptives.lookupSqrt(numAvailableInds);            // pre-calculate a square root lookup table
		Descriptives.initializeZScoreToPValue();
		
		boolean permuting = false;
		
		System.out.println("Will write output to dir: " + m_settings.outputReportsDir);
		
		int permStart = 0;
		int permEnd = m_settings.nrPermutationsFDR + 1;
		if (m_settings.startWithPermutation != null) {
			permStart = m_settings.startWithPermutation;
			permEnd = permStart + m_settings.nrPermutationsFDR + 1;
		}
		
		boolean hasResults = true;
		
		DoubleMatrixDataset<String, String>[] covariateData = null;
		if (dataHasCovariates) {
			covariateData = new DoubleMatrixDataset[m_gg.length];
			for (int d = 0; d < m_gg.length; d++) {
				covariateData[d] = m_gg[d].getCovariateData();
			}
		}
		
		for (int permutationRound = permStart; permutationRound < permEnd; permutationRound++) {
			RunTimer permtime = new RunTimer();
			
			if (permutationRound > 0) {
				System.out.print("Permuting data, round: " + permutationRound + " of " + m_settings.nrPermutationsFDR + "\n" + ConsoleGUIElems.LINE);
				
				for (int d = 0; d < m_gg.length; d++) {
//                    int[] indWGAOriginal = m_gg[d].getExpressionToGenotypeIdArray();
					m_gg[d].permuteSampleLables(m_settings.randomNumberGenerator);
					if (m_settings.permuteCovariates) {
						m_gg[d].permuteCovariates(m_settings.randomNumberGenerator);
					}

//                    int[] indWGAPerm = m_gg[d].getExpressionToGenotypeIdArray();
////                    int identical = 0;
////                    for (int i = 0; i < indWGAPerm.length; i++) {
////                        if (indWGAOriginal[i] == indWGAPerm[i]) {
////                            identical++;
////                        }
////                    }
////
//////                    System.out.println("After permuting: " + identical + " unchanged, " + indWGAOriginal.length + " total");
				}
				permuting = true;
			} else {
				System.out.print("Running real eQTL analysis\n" + ConsoleGUIElems.LINE);
			}
			
			int[][] expressionToGenotypeIds = new int[m_gg.length][0];
			for (int d = 0; d < m_gg.length; d++) {
				expressionToGenotypeIds[d] = m_gg[d].getExpressionToGenotypeIdArray();
			}
			
			LinkedBlockingQueue<WorkPackage> resultQueue = new LinkedBlockingQueue<WorkPackage>(250);
			ResultProcessorThread resultthread = new ResultProcessorThread(m_settings.nrThreads, resultQueue, m_settings.createBinaryOutputFiles,
					m_gg, m_settings, m_probeTranslationTable, permuting, permutationRound, m_snpList, m_probeList, m_workPackages);
			resultthread.setName("ResultProcessorThread");
			resultthread.start();
			
			// start production in advance
			LinkedBlockingQueue<WorkPackage> packageQueue = new LinkedBlockingQueue<WorkPackage>(100000);
			WorkPackageProducer producer = new WorkPackageProducer(packageQueue, m_workPackages, m_snpList, m_probeList, m_probeTranslationTable, m_snpTranslationTable, m_gg, snploaders, m_settings, permuting);
			producer.setName("WorkPackageProducerThread");
			producer.start();
			
			// run calculations
			for (int tnum = 0; tnum < pool.length; tnum++) {
				EQTLPlotter plotter = null;
				if (!permuting) {
					plotter = new EQTLPlotter(m_gg, m_settings, m_probeList, m_probeTranslationTable);
				}
				pool[tnum] = new CalculationThread(permutationRound, packageQueue, resultQueue, expressiondata, covariateData, m_probeTranslationTable, expressionToGenotypeIds, m_settings, plotter, m_settings.createBinaryOutputFiles, m_settings.useAbsoluteZScorePValue, m_settings.confineSNPsToSNPsPresentInAllDatasets);
				pool[tnum].setName("CalcThread-" + tnum);
				pool[tnum].start();
				
			}
			
			// kill the threads
			try {
				
				producer.join();

//                System.out.println("Joining calculation threads");
				for (int threadNum = 0; threadNum < m_settings.nrThreads; threadNum++) {
					pool[threadNum].join();
				}
				
				WorkPackage poison = new WorkPackage();
				poison.results = new Result(true);
				resultQueue.put(poison);
				resultthread.join();
				
			} catch (InterruptedException e) {
				System.err.println("Exception: Main Thread interrupted.");
			}
			System.out.print(ConsoleGUIElems.LINE);
			System.out.println("Round done. Elapsed time:\t" + permtime.getTimeDesc());
			System.out.println("");
			resultQueue.clear();
			packageQueue.clear();
			
			resultQueue = null;
			packageQueue = null;
			permtime = null;
			producer = null;
			resultthread = null;
			expressionToGenotypeIds = null;
			for (int i = 0; i < pool.length; i++) {
				pool[i] = null;
			}
			
			// check whether there were results..
			if (m_settings.createTEXTOutputFiles) {
				String fileName;
				if (permutationRound > 0) {
					fileName = m_settings.outputReportsDir + "PermutedEQTLsPermutationRound" + permutationRound + ".txt.gz";
				} else {
					fileName = m_settings.outputReportsDir + "eQTLs.txt.gz";
				}
				TextFile tf = new TextFile(fileName, TextFile.R);
				tf.readLine(); // skip header
				int lnCounter = 0;
				String line = tf.readLine();
				while (line != null) {
					lnCounter++;
					if (lnCounter > 1) {
						break;
					}
					line = tf.readLine();
				}
				tf.close();
				if (lnCounter == 0) {
					System.err.println("WARNING: QTL Mapping did not yield any results.");
					hasResults = false;
				}
			}
		}
		
		for (int d = 0; d < snploaders.length; d++) {
			snploaders[d].close();
		}
		
		
		if (!m_settings.skipFDRCalculation || (!m_settings.runOnlyPermutations && hasResults)) {
			if (m_settings.createTEXTOutputFiles && m_settings.nrPermutationsFDR > 0) {
				System.out.println("Calculating FDR:\n" + ConsoleGUIElems.LINE);
				FDR.calculateFDR(m_settings.outputReportsDir, m_settings.nrPermutationsFDR, m_settings.maxNrMostSignificantEQTLs, m_settings.fdrCutOff, m_settings.createQQPlot, null, null, m_settings.fdrType, m_settings.fullFdrOutput);
				
				if (m_settings.createDotPlot) {
					EQTLDotPlot edp = new EQTLDotPlot();
					try {
						edp.draw(m_settings.outputReportsDir + "/eQTLsFDR" + m_settings.fdrCutOff + ".txt", m_settings.outputReportsDir + "/DotPlot-FDR" + m_settings.fdrCutOff + ".pdf", EQTLDotPlot.Output.PDF); // "/eQTLsFDR" + fdrCutOff + ".txt", outputReportsDir + "/eQTLsFDR" + fdrCutOff + "DotPlot.png"
					} catch (DocumentException ex) {
						Logger.getLogger(MetaQTL3.class.getName()).log(Level.SEVERE, null, ex);
					}
					edp = null;
				}
				
			}
		} else {
			String reason = "";
			if (m_settings.skipFDRCalculation) {
				reason = "Defined in settings.";
			} else if (m_settings.runOnlyPermutations) {
				reason = "Only running permutations.";
			} else if (!hasResults) {
				reason = "No results for QTL mapping.";
			}
			System.out.println("Skipping FDR calculation. Reason: " + reason);
		}
		
		
		System.out.print(ConsoleGUIElems.DOUBLELINE);
		
		System.out.println("eQTL mapping elapsed:\t" + t.getTimeDesc() + "\n");
	}
	
	protected long determineSNPProbeCombinations() throws IOException {
		String loc = m_settings.outputReportsDir + "excludedSNPsBySNPProbeCombinationFilter.txt.gz";
		TextFile excludedSNPs = new TextFile(loc, TextFile.W);
		long maxNrTestsToPerform = 0;
		int[] midpoint = new int[m_probeList.length];
		byte[] chr = new byte[m_probeList.length];
		HashMap<Byte, ArrayList<Integer>> chrToProbe = new HashMap<Byte, ArrayList<Integer>>();
		System.out.println("- Calculating probe midpoint positions");
		HashSet<String> visitedProbes = new HashSet<String>();
		
		for (int p = 0; p < m_probeList.length; p++) {
			for (int d = 0; d < m_gg.length; d++) {
				if (m_probeTranslationTable.get(d, p) != -9 && !visitedProbes.contains(m_probeList[p])) {
					int pid = m_probeTranslationTable.get(d, p);
					int start = m_gg[d].getExpressionData().getChrStart()[pid];
					int stop = m_gg[d].getExpressionData().getChrStop()[pid];
					midpoint[p] = (int) Math.floor((double) (stop + start) / 2);
					chr[p] = m_gg[d].getExpressionData().getChr()[pid];
					
					ArrayList<Integer> probes = chrToProbe.get(chr[p]);
					if (probes == null) {
						probes = new ArrayList<Integer>();
					}
					probes.add(p);
					
					chrToProbe.put(chr[p], probes);
					visitedProbes.add(m_probeList[p]);
				}
			}
		}
		
		WorkPackage[] workPackages = new WorkPackage[m_snpList.length];
		// improve performance by sorting here, and then breaking later..
		int numWorkPackages = 0;
		
		System.out.println("- Determining SNP-Probe combinations to test");
		
		boolean cisOnly = false;
		boolean cisTrans = false;
		boolean transOnly = false;
		
		if (m_settings.cisAnalysis && !m_settings.transAnalysis) {
			cisOnly = true;
		} else if (!m_settings.cisAnalysis && m_settings.transAnalysis) {
			transOnly = true;
		} else if (m_settings.cisAnalysis && m_settings.transAnalysis) {
			cisTrans = true;
		}
		
		HashMap<String, Integer> probeNameToId = null;
		if (m_settings.tsSNPProbeCombinationsConfine != null) {
			m_settings.cisAnalysis = true;
			m_settings.transAnalysis = false;
			for (int i = 0; i < m_gg.length; i++) {
				m_gg[i].getSettings().cisAnalysis = true;
				m_gg[i].getSettings().transAnalysis = false;
			}
			cisTrans = false;
			cisOnly = true;
			transOnly = false;
			
			probeNameToId = new HashMap<String, Integer>();
			for (int i = 0; i < m_probeList.length; i++) {
				probeNameToId.put(m_probeList[i], i);
			}
		}
		
		int prevProc = 0;
		ProgressBar pb = new ProgressBar(m_snpList.length);
		for (int s = 0; s < m_snpList.length; s++) {
			WorkPackage output = null;
			
			SNP[] snps = new SNP[m_gg.length];
			
			byte snpchr = -1;
			int snppos = -1;
			int dreq = 0;
			
			String snpname = "";
			
			// load the genomic positions for this snp
			for (int d = 0; d < m_gg.length; d++) {
				Integer snpId = m_snpTranslationTable.getQuick(d, s);
				if (snpId != -9) {
					snps[d] = m_gg[d].getGenotypeData().getSNPObject(snpId);
					snpchr = snps[d].getChr();
					snppos = snps[d].getChrPos();
					snpname = snps[d].getName();
					dreq = d;
				}
			}
			
			ArrayList<Integer> probeOnChr = chrToProbe.get(snpchr);
			
			// cis trans
			if (cisTrans) {
				output = new WorkPackage();
				output.setMetaSNPId(s);
				output.setSnps(snps);
				output.setProbes(null);
				numWorkPackages++;
				maxNrTestsToPerform += m_probeList.length;
				workPackages[s] = output;
				
				// cis or trans
			} else {
				ArrayList<Integer> probeToTest = null;
				if (m_settings.tsSNPProbeCombinationsConfine != null) {
					HashSet<String> probesSelected;
					if (m_settings.snpProbeConfineBasedOnChrPos) {
						probesSelected = m_settings.tsSNPProbeCombinationsConfine.get(snpchr + ":" + snppos);
					} else {
						probesSelected = m_settings.tsSNPProbeCombinationsConfine.get(snpname);
					}
					
					if (probesSelected != null && probeNameToId != null) {
						probeToTest = new ArrayList<Integer>();
						for (String probe : probesSelected) {
							Integer probeId = probeNameToId.get(probe);
							if (probeId == null) {
								System.err.println("You selected the following SNP-Probe combination, but probe not present in dataset.\t" + snpname + "\t-\t" + probe);
							} else {
								probeToTest.add(probeId);
							}
						}
					}
				} else {
					if (probeOnChr != null && !probeOnChr.isEmpty()) {
						probeToTest = new ArrayList<Integer>();
						for (int e = 0; e < probeOnChr.size(); e++) {
							int p = probeOnChr.get(e);
							
							if (Math.abs(midpoint[p] - snppos) < m_settings.ciseQTLAnalysMaxSNPProbeMidPointDistance) {
								probeToTest.add(p); // depending if the test is cis or trans, we test these probes, or not.
								// System.out.println(snps[dreq].getName()+"\t"+p+"\t"+probeList[p]+"\t"+Math.abs(midpoint[p] - snppos));
							}
						}
					}
				}
				
				// don't add the cis-workpackage when there are no probes to test.
				// cisonly is also used when performing analysis on a certain set of snp-probe combos.
				if (cisOnly && (probeToTest == null || probeToTest.isEmpty())) {
					workPackages[s] = null;
					
					// reduce the output size of the snp filter output to query snps, if any
					if (m_settings.tsSNPsConfine == null || m_settings.tsSNPsConfine.contains(m_snpList[s])) {
						excludedSNPs.write(snpname + "\tNo probes to test.\n");
					}
				} else {
					int[] testprobes = new int[0];
					if (probeToTest != null) {
						testprobes = new int[probeToTest.size()];
						for (int p = 0; p < testprobes.length; p++) {
							testprobes[p] = probeToTest.get(p);
						}
					}
					
					output = new WorkPackage();
					output.setMetaSNPId(s);
					output.setSnps(snps);
					output.setProbes(testprobes);
					
					workPackages[s] = output;
					
					numWorkPackages++;
					if (cisOnly) {
						maxNrTestsToPerform += testprobes.length;
					} else {
						if (testprobes.length != 0) {
							maxNrTestsToPerform += (m_probeList.length - testprobes.length);
						} else {
							maxNrTestsToPerform += (m_probeList.length);
						}
						
					}
				}
			}
			
			pb.iterate();
		}
		
		pb.close();
		System.out.println("");
		if (numWorkPackages != m_snpList.length) {
			m_workPackages = new WorkPackage[numWorkPackages];
			int q = 0;
			for (int i = 0; i < workPackages.length; i++) {
				if (workPackages[i] != null) {
					m_workPackages[q] = workPackages[i];
					m_workPackages[q].setId(q);
					q++;
				}
			}
		} else {
			for (int i = 0; i < workPackages.length; i++) {
				workPackages[i].setId(i);
			}
			this.m_workPackages = workPackages;
		}

//		java.util.Arrays.sort(m_workPackages);
		
		excludedSNPs.close();
		
		System.out.println("The maximum number of SNPs to test: " + m_workPackages.length);
		System.out.println("The maximum number of SNP-Probe combinations: " + maxNrTestsToPerform);
		
		if (m_settings.maxNrMostSignificantEQTLs > maxNrTestsToPerform) {
			m_settings.maxNrMostSignificantEQTLs = (int) maxNrTestsToPerform;
		}
		
		return (maxNrTestsToPerform);
	}
	
	protected void printSummary() {
		
		System.out.print("\nSummary\n" + ConsoleGUIElems.DOUBLELINE);
		int totalSamples = 0;
		for (TriTyperGeneticalGenomicsDataset m_gg1 : m_gg) {
			System.out.print("Dataset:\t" + m_gg1.getSettings().name);
			System.out.print("\tprobes:\t" + m_gg1.getExpressionData().getProbes().length);
			System.out.print("\tSNPs:\t" + m_gg1.getGenotypeData().getSNPs().length);
			totalSamples += m_gg1.getTotalGGSamples();
			System.out.println("\tsamples:\t" + m_gg1.getTotalGGSamples());
		}
		System.out.println("");
		System.out.print("\nTotals\n" + ConsoleGUIElems.DOUBLELINE);
		System.out.println("Total number of datasets:\t" + m_gg.length);
		System.out.println("Total number of samples:\t" + totalSamples);
		System.out.println("Maximum number of SNPs to test:\t" + m_workPackages.length);
		System.out.println("Maximum number of Probes to test:\t" + m_probeList.length);
		
		if (totalSamples == 0) {
			System.err.println("ERROR!: No samples detected");
			System.exit(-1);
		}
		System.out.println("");
		System.out.print("\nAnalysis\n" + ConsoleGUIElems.DOUBLELINE);
		
		if (m_settings.cisAnalysis && m_settings.transAnalysis) {
			System.out.println("- cis/trans analysis");
		} else if (m_settings.cisAnalysis && !m_settings.transAnalysis) {
			System.out.println("- cis analysis");
		} else if (!m_settings.cisAnalysis && m_settings.transAnalysis) {
			System.out.println("- trans analysis");
		}
		if (m_settings.metaAnalyseInteractionTerms) {
			System.out.println("- interaction analysis");
			if (!m_settings.performParametricAnalysis) {
				System.out.println("- WARNING: running interaction model on non-parametric data!");
			}
		}
		if (!m_settings.performParametricAnalysis) {
			System.out.println("- non-parametric (Spearman ranked) correlation");
		} else {
			System.out.println("- parametric (Pearson) correlation");
		}
		System.out.println("- Mid-point distance:\t" + m_settings.ciseQTLAnalysMaxSNPProbeMidPointDistance);
		System.out.println("- FDR cutoff:\t" + m_settings.fdrCutOff);
		System.out.println("- Nr. permutations:\t" + m_settings.nrPermutationsFDR);
		System.out.println("- Nr. Threads:\t" + m_settings.nrThreads);
		System.out.println("- Max nr results:\t" + m_settings.maxNrMostSignificantEQTLs);
		
		if (m_settings.createBinaryOutputFiles) {
			System.out.println("- creating BINARY output");
		}
		
		if (m_settings.createTEXTOutputFiles) {
			System.out.println("- creating TEXT output");
		}
		
		System.out.println("");
	}
}
