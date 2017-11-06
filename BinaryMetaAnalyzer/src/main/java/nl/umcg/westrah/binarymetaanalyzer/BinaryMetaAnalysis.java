/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.westrah.binarymetaanalyzer;

import gnu.trove.map.hash.TObjectIntHashMap;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;

import static com.google.common.primitives.Ints.toArray;

/**
 * @author Harm-Jan
 */
public class BinaryMetaAnalysis {
	
	
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
			System.out.println("Usage of the binary meta-analysis: settings.xml replacetext replacetextwith");
			System.exit(-1);
		}
		
		BinaryMetaAnalysis meta = new BinaryMetaAnalysis(settingsFile, textToReplace, replaceTextWith);
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
	private boolean sorted;
	private int locationToStoreResult;
	private MetaQTL4MetaTrait[][] snpprobeCombos;
	
	TObjectIntHashMap<MetaQTL4MetaTrait> traitMap = new TObjectIntHashMap<MetaQTL4MetaTrait>();
	MetaQTL4MetaTrait[] traitList = null;
	
	public BinaryMetaAnalysis(String settingsFile, String textToReplace, String replaceTextWith) {
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
		
		finalEQTLs = new QTL[(maxResults + tmpbuffersize)];
		try {
			run();
		} catch (IOException ex) {
			Logger.getLogger(BinaryMetaAnalysis.class.getName()).log(Level.SEVERE, null, ex);
		}
		
	}
	
	public void run() throws IOException {
		
		String outdir = settings.getOutput();
		outdir = Gpio.formatAsDirectory(outdir);
		Gpio.createDir(outdir);
		// load probe annotation and index
		// this particular probe annotation can take multiple probes for a single location into account.
		System.out.println("Loading probe annotation from: " + settings.getProbetranslationfile());
		loadProbeAnnotation();
		
		if (traitList.length == 0) {
			System.err.println("Error: no annotation loaded.");
			System.exit(-1);
		}
		
		for (int permutation = settings.getStartPermutations(); permutation <= settings.getNrPermutations(); permutation++) {
			clearResultsBuffer();
			
			// create dataset objects
			System.out.println("Running permutation " + permutation);
			datasets = new BinaryMetaAnalysisDataset[settings.getDatasetlocations().size()];
			
			System.out.println("Loading datasets");
			for (int d = 0; d < datasets.length; d++) {
				datasets[d] = new BinaryMetaAnalysisDataset(settings.getDatasetlocations().get(d),
						settings.getDatasetnames().get(d),
						settings.getDatasetPrefix().get(d),
						permutation,
						settings.getDatasetannotations().get(d),
						probeAnnotation);
			}
			System.out.println("Loaded " + datasets.length + " datasets");
			
			// create meta-analysis SNP index. have to recreate this every permutation,
			// since the order of SNPs is generated at random.
			System.out.println("Creating SNP index");
			createSNPIndex(outdir);
			System.out.println("Total of " + snpIndex.length + " SNPs");
			System.out.println("Creating probe index");
			createProbeIndex(outdir);
			System.out.println("Total of " + probeIndex.length + " probes");
			
			// make index of snp/probe combinations, if any specified
			createSNPProbeCombos(outdir);
			
			// if snp probe selection is set, make combinations here..
			if (snpChr == null) {
				System.out.println("Loading SNP annotation from " + settings.getSNPAnnotationFile());
				loadSNPAnnotation();
			}
			
			// run analysis
			System.out.println("Type of analysis: " + settings.getAnalysisType());
			System.out.println("Cis-window: " + settings.getCisdistance());
			System.out.println("Trans-window: " + settings.getTransdistance());
			
			System.out.println("Starting meta-analysis");
			ProgressBar pb = new ProgressBar(snpList.length);
			for (int snp = 0; snp < snpList.length; snp++) {
				boolean printed = false;
				int[] sampleSizes = new int[datasets.length];
				Boolean[] flipZScores = new Boolean[datasets.length];
				String alleles = null;
				String alleleAssessed = null;
				
				// determine whether to flip the alleles for a certain dataset
				for (int d = 0; d < datasets.length; d++) {
					int datasetSNPId = snpIndex[snp][d];
					if (datasetSNPId != -9) {
						sampleSizes[d] = datasets[d].getSampleSize(datasetSNPId);
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
				
				// get ZScores for this SNP
				// get list of probes to test
				if (settings.getAnalysisType().equals(BinaryMetaAnalysisSettings.Analysis.CIS) || snpprobeCombos != null) {
					// do cis stuff, or stuff to specific sets of probes...
					
					HashMap<MetaQTL4MetaTrait, Integer> cisProbeMap = new HashMap<MetaQTL4MetaTrait, Integer>();
					MetaQTL4MetaTrait[] cisProbeArray = null;
					if (snpprobeCombos != null) {
						cisProbeArray = snpprobeCombos[snp];
					} else {
						Set<MetaQTL4MetaTrait> cisProbesForSNP = probeAnnotation.getMetatraits().getTraitInWindow(snpChr[snp], snpPositions[snp], settings.getCisdistance());
						cisProbeArray = cisProbesForSNP.toArray(new MetaQTL4MetaTrait[0]);
					}
					
					if (cisProbeArray == null || cisProbeArray.length == 0) {
						// nothing to do. skip variant //
					} else {
						int ctr = 0;
						for (MetaQTL4MetaTrait cisProbe : cisProbeArray) {
							cisProbeMap.put(cisProbe, ctr);
							ctr++;
						}
						
						
						float[][] finalZScores = new float[cisProbeMap.size()][datasets.length];
						
						// get list of probes to test for each dataset
						for (int d = 0; d < datasets.length; d++) {
							if (flipZScores[d] == null) {
								// the allele could not be flipped. set the Z to NaN
								for (float[] zScore : finalZScores) {
									zScore[d] = Float.NaN;
								}
							} else {
								//initialize z-score
								for (int p = 0; p < cisProbeMap.size(); p++) {
									finalZScores[p][d] = Float.NaN; // this is not very nice, but does prevent the metaZ method from going nuts
								}
								
								// load the z-scores for the dataset
								int datasetSNPId = snpIndex[snp][d];
								
								if (datasetSNPId != -9) { // -9 means: snp not available
									float[] datasetZScores = datasets[d].getZScores(datasetSNPId);
									
									if (datasets[d].getIsCisDataset()) {
										// this requires us to retrieve the z-scores differently:
										// a cis dataset only stores z-scores for the tested probes/traits
										// position 0 in the datasetZScores array may therefore point to meta-trait 1000 in our annotation
										// we need to figure out which probes match up.. their orders might be different
										// and the number of probes tested in each dataset might differ as well
										
										// get the probes tested against the SNP
										MetaQTL4MetaTrait[] datasetCisProbes = datasets[d].getCisProbes(datasetSNPId);
										
										for (int i = 0; i < datasetCisProbes.length; i++) {
											MetaQTL4MetaTrait p = datasetCisProbes[i];
											if (p != null) {
												Integer index = cisProbeMap.get(p);
												if (index != null) {
													float datasetZ = datasetZScores[i];
													finalZScores[index][d] = datasetZ;
													if (flipZScores[d]) {
														finalZScores[index][d] *= -1;
													}
												}
											}
										}
									} else { // this is not a cis dataset
										// use the full probe index
										for (int probe = 0; probe < cisProbeArray.length; probe++) {
											MetaQTL4MetaTrait cisProbe = cisProbeArray[probe];
											Integer metaProbeIndex = traitMap.get(cisProbe);
											Integer datasetProbeId = probeIndex[metaProbeIndex][d];
											if (datasetProbeId != null) {
												finalZScores[probe][d] = datasetZScores[datasetProbeId];
												if (flipZScores[d]) {
													finalZScores[probe][d] *= -1;
												}
											}
										}
									}
									
									
								}
							}
						}
						
						// meta-analyze!
						for (int probe = 0; probe < finalZScores.length; probe++) {
							MetaQTL4MetaTrait t = cisProbeArray[probe];
							
							double metaZ = ZScores.getWeightedZ(finalZScores[probe], sampleSizes);
							double p = Descriptives.convertZscoreToPvalue(metaZ);
							
							if (!Double.isNaN(p) && !Double.isNaN(metaZ)) {
								// create output object
								QTL q = new QTL(p, t, snp, BaseAnnot.toByte(alleleAssessed), metaZ, BaseAnnot.toByteArray(alleles), finalZScores[probe], sampleSizes); // sort buffer if needed.
//                            System.out.println(q.getSNPId()+"\t"+q.getMetaTrait().getMetaTraitName()+"\t"+q.toString());
								addEQTL(q);
							} else {
//                            if (!printed) {
//                                printed = true;
//                                System.out.println("SNP creates NaN pval: " + snpList[snp] + "\n");//z " + Strings.concat(zScores, Strings.semicolon) + "\tn " + Strings.concat(sampleSizes, Strings.semicolon));
//                            }
							}
						}
					}
				} else {
					// analysis is not cis, but may be cis/trans
					Set<MetaQTL4MetaTrait> cisProbes = null;
					if (!settings.getAnalysisType().equals(BinaryMetaAnalysisSettings.Analysis.CISTRANS)) {
						// do not test the cis probes if not cistrans
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
					for (int probe = 0; probe < traitList.length; probe++) {
						MetaQTL4MetaTrait t = traitList[probe];
						double metaAnalysisZ = ZScores.getWeightedZ(finalZScores[probe], sampleSizes);
						double metaAnalysisP = Descriptives.convertZscoreToPvalue(metaAnalysisZ);
						
						// create output object
						if (!Double.isNaN(metaAnalysisP) && !Double.isNaN(metaAnalysisZ)) {
							QTL q = new QTL(metaAnalysisP, t, snp, BaseAnnot.toByte(alleleAssessed), metaAnalysisZ, BaseAnnot.toByteArray(alleles), finalZScores[probe], sampleSizes); // sort buffer if needed.
//                            System.out.println(q.getSNPId()+"\t"+q.getMetaTrait().getMetaTraitName()+"\t"+q.toString());
							addEQTL(q);
						}
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
	
	private void createSNPProbeCombos(String outdir) throws IOException {
		
		snpprobeCombos = null;
		if (settings.getSNPProbeSelection() != null) {
			System.out.println("Getting SNPs from SNP/Probe selection file: " + settings.getSNPProbeSelection());
			TextFile tf = new TextFile(settings.getSNPProbeSelection(), TextFile.R);
			
			
			// make a quick snpmap
			HashMap<String, Integer> snpMap = new HashMap<String, Integer>();
			for (int i = 0; i < snpList.length; i++) {
				snpMap.put(snpList[i], i);
			}
			
			// make a quick traitmap
			HashMap<String, Integer> tmpTraitMap = new HashMap<String, Integer>();
			for (int i = 0; i < traitList.length; i++) {
				tmpTraitMap.put(traitList[i].getMetaTraitName(), i);
			}
			
			// combine
			ArrayList<SNPProbePair> pairs = new ArrayList<SNPProbePair>();
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				if (elems.length >= 2) {
					String snp = elems[0];
					String trait = elems[1];
					
					// find the snp
					Integer snpId = snpMap.get(snp);
					if (snpId != null) {
						// try to find the probe
						Integer traitId = tmpTraitMap.get(trait);
						if (traitId != null) {
							SNPProbePair p = new SNPProbePair(snpId, traitId);
							pairs.add(p);
						}
					}
				}
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
			
			if (pairs.isEmpty()) {
				System.err.println("ERROR: SNP/Probe file defined, but none of the specified combinations found in data.");
				System.exit(-1);
			}
			
			
			Collections.sort(pairs);
			
			snpprobeCombos = new MetaQTL4MetaTrait[snpList.length][];
			
			int prevSNP = -1;
			ArrayList<MetaQTL4MetaTrait> list = new ArrayList<MetaQTL4MetaTrait>();
			for (SNPProbePair p : pairs) {
				if (p.getSnpId() != prevSNP) {
					if (prevSNP > -1) {
						snpprobeCombos[prevSNP] = list.toArray(new MetaQTL4MetaTrait[0]);
						list = new ArrayList<MetaQTL4MetaTrait>();
					}
				}
				list.add(traitList[p.getProbeId()]);
				prevSNP = p.getSnpId();
			}
		}
	}
	
	private void createSNPIndex(String outdir) throws IOException {
		
		HashSet<String> confineToTheseSNPs = null;
		
		HashSet<String> snpPreSelection = null;
		if (settings.getSNPProbeSelection() != null) {
			System.out.println("Getting SNPs from SNP/Probe selection file: " + settings.getSNPProbeSelection());
			snpPreSelection = new HashSet<String>();
			TextFile tf = new TextFile(settings.getSNPProbeSelection(), TextFile.R);
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				String snp = elems[0];
				snpPreSelection.add(snp);
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
			System.out.println("Found " + snpPreSelection.size() + " unique snps in SNP/Probe selection file.");
			if (snpPreSelection.isEmpty()) {
				System.err.println("Error: SNP/Probe selection file defined, but no SNPs found.");
				System.exit(-1);
			}
		}
		
		if (settings.getSNPSelection() != null) {
			System.out.println("Selecting SNPs from file: " + settings.getSNPSelection());
			confineToTheseSNPs = new HashSet<String>();
			TextFile tf = new TextFile(settings.getSNPSelection(), TextFile.R);
			ArrayList<String> snps = tf.readAsArrayList();
			tf.close();
			if (snpPreSelection == null) {
				confineToTheseSNPs.addAll(snps);
			} else {
				System.out.println("Intersecting with SNP/Probe selection.");
				for (String snp : snps) {
					if (snpPreSelection.contains(snp)) {
						confineToTheseSNPs.add(snp);
					}
				}
			}
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
		
		TextFile tf = new TextFile(outdir + "snpindex.txt", TextFile.W);
		String header = "metaID";
		for (int d = 0; d < datasets.length; d++) {
			header += "\t" + datasets[d].getName() + "-sid";
		}
		tf.writeln(header);
		for (int s = 0; s < snpList.length; s++) {
			String ln = snpList[s];
			for (int d = 0; d < datasets.length; d++) {
				ln += "\t" + snpIndex[s][d];
			}
			tf.writeln(ln);
		}
		tf.close();
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
	private void createProbeIndex(String outdir) throws IOException {
		
		HashSet<String> confineToTheseProbes = null;
		
		HashSet<String> probePreselection = null;
		if (settings.getSNPProbeSelection() != null) {
			if (settings.getSNPProbeSelection() != null) {
				System.out.println("Getting Probes from SNP/Probe selection file: " + settings.getSNPProbeSelection());
				probePreselection = new HashSet<String>();
				TextFile tf = new TextFile(settings.getSNPProbeSelection(), TextFile.R);
				String[] elems = tf.readLineElems(TextFile.tab);
				while (elems != null) {
					if (elems.length >= 2) {
						String probe = elems[1];
						probePreselection.add(probe);
					}
					elems = tf.readLineElems(TextFile.tab);
				}
				tf.close();
				System.out.println("Found " + probePreselection.size() + " unique probes in SNP/Probe selection file.");
				if (probePreselection.isEmpty()) {
					System.err.println("Error: SNP/Probe selection file defined, but no Probes found.");
					System.exit(-1);
				}
			}
		}
		
		if (settings.getProbeselection() != null) {
			System.out.println("Selecting Probes from file: " + settings.getProbeselection());
			confineToTheseProbes = new HashSet<String>();
			TextFile tf = new TextFile(settings.getProbeselection(), TextFile.R);
			
			if (probePreselection == null) {
				confineToTheseProbes.addAll(tf.readAsArrayList());
			} else {
				ArrayList<String> confineTMP = tf.readAsArrayList();
				for (String p : confineTMP) {
					if (probePreselection.contains(p)) {
						confineToTheseProbes.add(p);
					}
				}
			}
			tf.close();
			System.out.println(confineToTheseProbes.size() + " Probes loaded.");
		}
		
		
		System.out.println("");
		probeIndex = new Integer[traitList.length][datasets.length];
		
		for (int d = 0; d < datasets.length; d++) {
			String[] probes = datasets[d].getProbeList();
			int platformId = probeAnnotation.getPlatformId(datasets[d].getPlatform());
			
			HashMap<String, MetaQTL4MetaTrait> traitHashForPlatform = probeAnnotation.getTraitHashForPlatform(platformId);
			System.out.println(probeAnnotation.getTraitHashPerPlatform().size());
			
			System.out.println(datasets[d].getName() + "\t" + platformId + "\t" + datasets[d].getPlatform() + "\t" + traitHashForPlatform.size());
			for (int p = 0; p < probes.length; p++) {
				
				MetaQTL4MetaTrait t = traitHashForPlatform.get(probes[p]);
				int index = traitMap.get(t);
				
				if (confineToTheseProbes == null || confineToTheseProbes.contains(probes[p]) || confineToTheseProbes.contains(t.getMetaTraitName())) {
					// TODO: selecting on individual platform probes/genes. Is this the best approach?
					probeIndex[index][d] = p;
				}
			}
		}
		
		System.out.println("");
		
		TextFile out = new TextFile(outdir + "probeindex.txt", TextFile.W);
		
		String header = "metaID";
		for (int d = 0; d < datasets.length; d++) {
			header += "\t" + datasets[d].getName() + "-pid\t" + datasets[d].getName() + "-probename";
		}
		out.writeln(header);
		for (int p = 0; p < probeIndex.length; p++) {
			
			String lnout = "" + traitList[p].getMetaTraitId();
			for (int d = 0; d < datasets.length; d++) {
				Integer pid = probeIndex[p][d];
				String probeName = null;
				if (pid != null) {
					probeName = datasets[d].getProbeList()[pid];
				}
				lnout += "\t" + pid + "\t" + probeName;
			}
			
			out.writeln(lnout);
		}
		
		out.close();
	}
	
	private void addEQTL(QTL q) {
		
		double pval = q.getPvalue();
		if (bufferHasOverFlown) {
			if (pval <= maxSavedPvalue) {
				
				sorted = false;
				
				finalEQTLs[locationToStoreResult] = q;
				locationToStoreResult++;
				
				if (locationToStoreResult == finalEQTLs.length) {
					
					Arrays.sort(finalEQTLs);
					sorted = true;
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
				String[] dsZBuilder = new String[datasets.length];
				
				for (int d = 0; d < datasetZScores.length; d++) {
					
					if (!Float.isNaN(datasetZScores[d])) {
						String str = format.format(datasetZScores[d]);
						
						dsBuilder[d] = settings.getDatasetnames().get(d);
						dsNBuilder[d] = "" + q.getDatasetSampleSizes()[d];
						dsZBuilder[d] = str;
					} else {
						dsBuilder[d] = "-";
						dsNBuilder[d] = "-";
						dsZBuilder[d] = "-";
					}
				}
				
				sb.append("\t");
				sb.append(Strings.concat(dsBuilder, Strings.semicolon));
				
				sb.append("\t");
				sb.append(Strings.concat(dsZBuilder, Strings.semicolon));
				
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
	
	private class SNPProbePair implements Comparable<SNPProbePair> {
		private int snpId;
		private int probeId;
		
		public SNPProbePair(int snpId, int probeId) {
			this.snpId = snpId;
			this.probeId = probeId;
		}
		
		public int getSnpId() {
			return snpId;
		}
		
		public void setSnpId(int snpId) {
			this.snpId = snpId;
		}
		
		public int getProbeId() {
			return probeId;
		}
		
		public void setProbeId(int probeId) {
			this.probeId = probeId;
		}
		
		@Override
		public boolean equals(Object o) {
			if (this == o) return true;
			if (o == null || getClass() != o.getClass()) return false;
			
			SNPProbePair that = (SNPProbePair) o;
			
			if (snpId != that.snpId) return false;
			return probeId == that.probeId;
		}
		
		@Override
		public int hashCode() {
			int result = snpId;
			result = 31 * result + probeId;
			return result;
		}
		
		
		@Override
		public int compareTo(SNPProbePair o) {
			if (this.equals(o)) {
				return 0;
			} else if (this.snpId == o.snpId) {
				if (this.probeId > o.probeId) {
					return 1;
				} else {
					return -1;
				}
			} else if (this.snpId < o.snpId) {
				return -1;
			} else {
				return 1;
			}
		}
	}
}
