/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.westrah.binarymetaanalyzer;

import gnu.trove.map.hash.TObjectIntHashMap;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * @author Harm-Jan
 */
public class InternalMetaAnalysis {

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
			System.out.println("Usage of the internal meta-analysis: settings.xml");
			
			System.exit(-1);
		}

		InternalMetaAnalysis meta = new InternalMetaAnalysis(settingsFile, textToReplace, replaceTextWith);
		System.exit(0);

	}

	private MetaQTL4TraitAnnotation probeAnnotation;
	private BinaryMetaAnalysisDataset dataset;
	private final InternalMetaAnalysisSettings settings;
	
	TObjectIntHashMap<MetaQTL4MetaTrait> traitMap = new TObjectIntHashMap<MetaQTL4MetaTrait>();
	MetaQTL4MetaTrait[] traitList = null;

	public InternalMetaAnalysis(String settingsFile, String textToReplace, String replaceTextWith) {
		// initialize settings
		settings = new InternalMetaAnalysisSettings();
		settings.parse(settingsFile, textToReplace, replaceTextWith);

		try {
			run();
		} catch (IOException ex) {
			Logger.getLogger(InternalMetaAnalysis.class.getName()).log(Level.SEVERE, null, ex);
		}

	}

	public void run() throws IOException {

		String outdir = settings.getOutput();

		for (int permutation = settings.getStartPermutations(); permutation <= settings.getNrPermutations(); permutation++) {
			clearResultsBuffer();

			// create dataset objects
			System.out.println("Running permutation " + permutation);
			dataset =  new BinaryMetaAnalysisDataset(settings.getDatasetlocation(),
						settings.getDatasetname(),
						settings.getDatasetPrefix(),
						permutation);

			// create meta-analysis SNP index. have to recreate this every permutation,
			// since the order of SNPs is generated at random.
			System.out.println("Creating SNP index");
			createSNPIndex(outdir);
			System.out.println("Total of " + snpIndex.length + " SNPs");
			System.out.println("Creating probe index");
			createProbeIndex(outdir);
			System.out.println("Total of " + probeIndex.length + " probes");

			if (snpChr == null) {
				System.out.println("Loading SNP annotation from " + settings.getSNPAnnotationFile());
				loadSNPAnnotation();
			}

			// run analysis
			System.out.println("Cis-analysis: " + settings.isCis());
			System.out.println("Trans-analysis: " + settings.isTrans());
			System.out.println("Cis-window: " + settings.getCisdistance());
			System.out.println("Trans-window: " + settings.getTransdistance());

			TextFile zscoreTableTf = null;

			if (settings.isMakezscoretable()) {

				String tableoutfile = outdir + "ZScoreMatrix-Permutation" + permutation + ".txt.gz";
				if (permutation == 0) {
					tableoutfile = outdir + "ZScoreMatrix.txt.gz";
				}
				System.out.println("Writing z-score table: " + tableoutfile);
				zscoreTableTf = new TextFile(tableoutfile, TextFile.W);

				// write header

				zscoreTableTf.writeln(zscoretableheader);
			}

			System.out.println("Starting meta-analysis");
			ProgressBar pb = new ProgressBar(snpList.length);
			for (int snp = 0; snp < snpList.length; snp++) {

				double[] zscoretableoutput = null;
				if (settings.isMakezscoretable()) {

					if (zscoretableoutput != null) {
						for (int i = 0; i < zscoretableoutput.length; i++) {
							zscoretableoutput[i] = Double.NaN;
						}
					} else {
						zscoretableoutput = new double[probeIndex.length];
					}
				}

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

				// get ZScores for this SNP, no matter what.
				// get list of probes to test
				if (settings.isCis() && !settings.isTrans()) {
					// do cis stuff

					// get all the possible traits near the SNP
					Set<MetaQTL4MetaTrait> cisProbesForSNP = probeAnnotation.getMetatraits().getTraitInWindow(snpChr[snp], snpPositions[snp], settings.getCisdistance());
					MetaQTL4MetaTrait[] cisProbeArray = cisProbesForSNP.toArray(new MetaQTL4MetaTrait[0]);
					HashMap<MetaQTL4MetaTrait, Integer> cisProbeMap = new HashMap<MetaQTL4MetaTrait, Integer>();
					int ctr = 0;
					for (MetaQTL4MetaTrait cisProbe : cisProbesForSNP) {
						cisProbeMap.put(cisProbe, ctr);
						ctr++;
					}

					float[][] finalZScores = new float[cisProbeMap.size()][datasets.length]; // size: [possible cis-probes][nr of datasets]

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
// this requires us to retrieve the z-scores differently
									// we need to figure out which probes match up, but their orders might be different
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


					for (int probe = 0; probe < finalZScores.length; probe++) {
						MetaQTL4MetaTrait t = cisProbeArray[probe];

						double metaZ = ZScores.getWeightedZ(finalZScores[probe], sampleSizes);
						double p = Descriptives.convertZscoreToPvalue(metaZ);

						if (settings.isMakezscoretable()) {
							// get the correct index for trait t
							int metaid = t.getCurrentMetaId();
							zscoretableoutput[metaid] = metaZ;
						}

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
				} else {
					Set<MetaQTL4MetaTrait> cisProbes = null;
					if (!settings.isCis()) {
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
					if (settings.isMakezscoretable()) {
						// initialize with NaN
						for (int i = 0; i < zscoretableoutput.length; i++) {
							zscoretableoutput[i] = Double.NaN;
						}
					}
					for (int probe = 0; probe < traitList.length; probe++) {
						MetaQTL4MetaTrait t = traitList[probe];
						double metaAnalysisZ = ZScores.getWeightedZ(finalZScores[probe], sampleSizes);
						double metaAnalysisP = Descriptives.convertZscoreToPvalue(metaAnalysisZ);

						if (settings.isMakezscoretable()) {
							zscoretableoutput[probe] = metaAnalysisZ;
						}

						// create output object
						if (!Double.isNaN(metaAnalysisP) && !Double.isNaN(metaAnalysisZ)) {
							QTL q = new QTL(metaAnalysisP, t, snp, BaseAnnot.toByte(alleleAssessed), metaAnalysisZ, BaseAnnot.toByteArray(alleles), finalZScores[probe], sampleSizes); // sort buffer if needed.
//                            System.out.println(q.getSNPId()+"\t"+q.getMetaTrait().getMetaTraitName()+"\t"+q.toString());
							addEQTL(q);
						}
					}
				}
				pb.iterate();

				// write z-score output
				if (settings.isMakezscoretable()) {
					String snpName = snpList[snp];
					// get alleles

					zscoreTableTf.writeln(snpName + "\t" + alleles + "\t" + alleleAssessed + "\t" + Strings.concat(zscoretableoutput, Strings.tab));
				}
			}
			pb.close();

			if (settings.isMakezscoretable()) {
				zscoreTableTf.close();
			}

			for (BinaryMetaAnalysisDataset dataset : datasets) {
				dataset.close();
			}
			writeBuffer(outdir, permutation);
		}
	}

	private void createSNPIndex(String outdir) throws IOException {

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

	// index the probes
	private void createProbeIndex(String outdir) throws IOException {

		HashSet<String> confineToTheseProbes = null;
		if (settings.getProbeselection() != null) {
			System.out.println("Selecting Probes from file: " + settings.getProbeselection());
			confineToTheseProbes = new HashSet<String>();
			TextFile tf = new TextFile(settings.getProbeselection(), TextFile.R);
			confineToTheseProbes.addAll(tf.readAsArrayList());
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

				if (probes[p].equals("60437")) {
					if (t != null) {
						System.out.println(t.getMetaTraitId());
					} else {
						System.out.println("not found");
					}
				}

				if (confineToTheseProbes == null || confineToTheseProbes.contains(probes[p])) {
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
        
        
        
//        public void run() {
//        nrProcessed = 0;
//        try {
//            if (m_createBinaryFiles) {
//                zScoreBinaryFile = new BinaryFile[m_gg.length];
//                zScoreRowNamesFile = new TextFile[m_gg.length];
//                if (m_gg.length > 1) {
//                    String metaAnalysisFileName = m_outputdir + "MetaAnalysis";
//                    if (m_permuting) {
//                        metaAnalysisFileName += "-PermutationRound-" + m_permutationround;
//                    }
//                    zScoreMetaAnalysisFile = new BinaryFile(metaAnalysisFileName + ".dat", BinaryFile.W);
//                    // write magic number
//                    if (m_cisOnly) {
//                        zScoreMetaAnalysisFile.writeInt(1);
//                    } else {
//                        zScoreMetaAnalysisFile.writeInt(0);
//                    }
//
//                    zScoreMetaAnalysisRowNamesFile = new TextFile(metaAnalysisFileName + "-RowNames.txt.gz", TextFile.W);
//                    zScoreMetaAnalysisRowNamesFile.writeln("SNP\tAlleles\tMinorAllele\tAlleleAssessed\tNrCalled");
//                    TextFile tf = new TextFile(metaAnalysisFileName + "-ColNames.txt.gz", TextFile.W);
//                    tf.writeList(Arrays.asList(m_probeList));
//                    tf.close();
//                }
//                for (int d = 0; d < m_gg.length; d++) {
//                    String fileName = m_outputdir + m_gg[d].getSettings().name;
//                    if (m_permuting) {
//                        fileName += "-PermutationRound-" + m_permutationround;
//                    }
//                    zScoreBinaryFile[d] = new BinaryFile(fileName + ".dat", BinaryFile.W);
//                    // write magic number
//                    if (m_cisOnly) {
//                        zScoreBinaryFile[d].writeInt(1);
//                    } else {
//                        zScoreBinaryFile[d].writeInt(0);
//                    }
//
//                    TextFile tf = new TextFile(fileName + "-ColNames.txt.gz", TextFile.W);
//                    tf.writeList(Arrays.asList(m_probeList));
//                    tf.close();
//                    zScoreRowNamesFile[d] = new TextFile(fileName + "-RowNames.txt.gz", TextFile.W);
//                    zScoreRowNamesFile[d].writeln("SNP\tAlleles\tMinorAllele\tAlleleAssessed\tNrCalled\tMaf\tHWE\tCallRate");
//                }
//            }
//
//            ProgressBar progressbar = new ProgressBar(m_availableWorkPackages.length);
//            boolean poison = false;
//
//            while (!poison) {
//                WorkPackage wp = m_queue.take();
//                Result r = wp.results;
//                if (wp.getHasResults()) {
//                    nrSNPsTested++;
//                }
//
//                if (r.poison) {
//                    poison = true;
//                } else if (r.pvalues != null) {
//
//                    nrTestsPerformed += wp.getNumTested();
//
//                    double[] pvalues = r.pvalues;
//
//                    //Is this working?
//                    if (m_createBinaryFiles && !poison) {
//                        writeBinaryResult(r);
//                    }
//
//                    if (m_createTEXTFiles && !poison) {
//                        // classic textual output.
//
//                        for (int p = 0; p < pvalues.length; p++) {
//                            double pval = pvalues[p];
//
//                            if (!Double.isNaN(pval) && pval <= highestP) {
//                                double[][] corr = r.correlations;
//                                double[] correlations = new double[corr.length];
//                                double[] zscores = new double[corr.length];
//                                int[] samples = new int[corr.length];
//
//                                double[] fc = new double[corr.length];
//                                double[] beta = new double[corr.length];
//                                double[] betase = new double[corr.length];
//
//                                for (int d = 0; d < correlations.length; d++) {
//                                    if (Double.isNaN(corr[d][p])) {
//                                        correlations[d] = Double.NaN;
//                                        zscores[d] = Double.NaN;
//                                        samples[d] = -9;
//                                        fc[d] = Double.NaN;
//                                        beta[d] = Double.NaN;
//                                        betase[d] = Double.NaN;
//                                    } else {
//                                        correlations[d] = corr[d][p];
//                                        if (m_useAbsoluteZScore) {
//                                            zscores[d] = Math.abs(r.zscores[d][p]);
//                                        } else {
//                                            zscores[d] = r.zscores[d][p];
//                                        }
//
//                                        samples[d] = r.numSamples[d];
//                                        fc[d] = r.fc[d][p];
//                                        beta[d] = r.beta[d][p];
//                                        betase[d] = r.se[d][p];
//                                    }
//                                }
////
//                                byte allele = -1;
//                                byte[] alleles = null;
//                                SNP[] snps = wp.getSnps();
//                                for (int d = 0; d < snps.length; d++) {
//                                    if (snps[d] != null) {
//                                        allele = snps[d].getMinorAllele();
//                                        alleles = snps[d].getAlleles();
//                                        break;
//                                    }
//                                }
//
//                                if (alleles == null) {
//                                    System.err.println("SNP has null alleles: ");
//                                    for (int d = 0; d < snps.length; d++) {
//
//                                        if (snps[d] != null) {
//
//                                            allele = snps[d].getMinorAllele();
//                                            System.err.println(allele);
//                                            alleles = snps[d].getAlleles();
//                                            System.err.println(alleles);
//                                            break;
//                                        }
//                                    }
//                                }
//
//                                double Zfinal = r.finalZScore[p];
//                                double finalbeta = r.finalBeta[p];
//                                double finalbetase = r.finalBetaSe[p];
//                                int pid;
//                                if (m_cisOnly) {
//                                    pid = wp.getProbes()[p];
//                                } else {
//                                    pid = p;
//                                }
//
//                                addEQTL(pid, wp.getId(), pval, Zfinal, correlations, zscores, samples, alleles, allele, fc, beta, betase, finalbeta, finalbetase);
//
//                            }
//                        }
//                    }
//
//                }
//
//                if (wp.results != null) {
//                    wp.clearResults();
//
//                }
//
//                progressbar.iterate();
//            }
//
//            progressbar.close();
//
//            //Is this working?
//            if (m_createBinaryFiles) {
//
//                String fileName = "check";
//                if (m_permuting) {
//                    fileName += "-PermutationRound-" + m_permutationround;
//                }
//                fileName += ".md5";
//
//                HexBinaryAdapter md5Parser = new HexBinaryAdapter();
//
//                BufferedWriter md5writer = new BufferedWriter(new FileWriter(m_outputdir + fileName));
//
//                for (int d = 0; d < m_gg.length; d++) {
//                    zScoreBinaryFile[d].close();
//
//                    fileName = m_gg[d].getSettings().name;
//                    if (m_permuting) {
//                        fileName += "-PermutationRound-" + m_permutationround;
//                    }
//                    fileName += ".dat";
//                    md5writer.write(md5Parser.marshal(zScoreBinaryFile[d].getWrittenHash()) + "  " + fileName + '\n');
//
//                    zScoreRowNamesFile[d].close();
//                }
//                if (m_gg.length > 1) {
//                    zScoreMetaAnalysisFile.close();
//
//                    fileName = "MetaAnalysis";
//                    if (m_permuting) {
//                        fileName += "-PermutationRound-" + m_permutationround;
//                    }
//                    fileName += ".dat";
//                    md5writer.write(md5Parser.marshal(zScoreMetaAnalysisFile.getWrittenHash()) + "  " + fileName + '\n');
//
//                    zScoreMetaAnalysisRowNamesFile.close();
//                }
//
//                md5writer.close();
//
//            }
//
//            if (m_createTEXTFiles) {
//                if (!sorted) {
//                    if (locationToStoreResult != 0) {
//
//                        Arrays.parallelSort(finalEQTLs, 0, locationToStoreResult);
////                        SmoothSort.sort(finalEQTLs, 0, locationToStoreResult);
////                        inplaceArrayQuickSort.sort(finalEQTLs, 0, locationToStoreResult);
//
//                    }
//                }
//                writeTextResults();
//            }
//
//        } catch (IOException e1) {
//            e1.printStackTrace();
//        } catch (InterruptedException e2) {
//            e2.printStackTrace();
//        }
//    }
//
//    private void writeBinaryResult(Result r) throws IOException {
//
//        if (r != null) {
//            int[] numSamples = null;
//            try {
//                numSamples = r.numSamples;
//            } catch (NullPointerException e) {
//                System.out.println("ERROR: null result?");
//            }
//
//            int wpId = r.wpid;
//            WorkPackage currentWP = m_availableWorkPackages[wpId];
//            double[][] zscores = r.zscores;
//
//            if (zscores != null) {
//                SNP[] snps = currentWP.getSnps();
//                int numDatasets = zscores.length;
//                double[] finalZscores = r.finalZScore;
//                StringBuilder snpoutput = null;
//
//                // if we're doing a meta-analysis, write the meta-analysis Z to a separate binaryFile
//                if (m_gg.length > 1) {
//                    int totalSampleNr = 0;
//                    String snpname = null;
//                    for (int d = 0; d < numDatasets; d++) {
//                        if (snps[d] != null) {
//                            snpname = snps[d].getName();
//
//                            byte[] alleles = snps[d].getAlleles();
//                            byte minorAllele = snps[d].getMinorAllele();
//                            byte alleleassessed = alleles[1];
//
//                            if (currentWP.getFlipSNPAlleles()[d]) {
//                                alleleassessed = alleles[0];
//                            }
//                            if (snpoutput == null) {
//                                snpoutput = new StringBuilder();
//                                snpoutput.append(snpname);
//                                snpoutput.append("\t");
//                                snpoutput.append(BaseAnnot.getAllelesDescription(alleles));
//                                snpoutput.append("\t");
//                                snpoutput.append(BaseAnnot.toString(minorAllele));
//                                snpoutput.append("\t");
//                                snpoutput.append(BaseAnnot.toString(alleleassessed));
//                            }
//                            totalSampleNr += r.numSamples[d];
//                        }
//                    }
//
//                    StringBuilder sb = null;
//                    for (int p = 0; p < finalZscores.length; p++) {
//                        float z = (float) finalZscores[p];
//                        if (m_cisOnly) {
//                            int[] probes = currentWP.getProbes();
//                            int probeId = probes[p];
//                            String probeName = m_probeList[probeId];
//                            if (sb == null) {
//                                sb = new StringBuilder();
//                            } else {
//                                sb.append("\t");
//                            }
//                            sb.append(probeName);
//
//                            zScoreMetaAnalysisFile.writeFloat(z);
//                        } else {
//                            zScoreMetaAnalysisFile.writeFloat(z);
//                        }
//                    }
//
//                    if (snpoutput != null) {
//                        snpoutput.append("\t");
//                        snpoutput.append(totalSampleNr);
//                        snpoutput.append("\t-\t-\t-\t");
//                        snpoutput.append(finalZscores.length);
//                        snpoutput.append("\t");
//                        if (sb != null) {
//                            snpoutput.append(sb.toString());
//                        } else {
//                            snpoutput.append("-");
//                        }
//                        zScoreMetaAnalysisRowNamesFile.writeln(snpoutput.toString());
//                    }
//                }
//
//                for (int d = 0; d < numDatasets; d++) {
//                    double[] datasetZScores = zscores[d];
//                    SNP datasetSNP = snps[d];
//                    if (datasetSNP != null) {
//                        BinaryFile outfile = zScoreBinaryFile[d];
//
//                        String snpname = datasetSNP.getName();
//
//                        byte[] alleles = datasetSNP.getAlleles();
//                        byte minorAllele = datasetSNP.getMinorAllele();
//                        byte alleleassessed = alleles[1];
//                        double hwe = datasetSNP.getHWEP();
//                        double cr = datasetSNP.getCR();
//                        double maf = datasetSNP.getMAF();
//
//                        if (currentWP.getFlipSNPAlleles()[d]) {
//                            alleleassessed = alleles[0];
//                        }
//                        TextFile snpfile = zScoreRowNamesFile[d];
//                        StringBuilder sb = null;
//                        for (int p = 0; p < datasetZScores.length; p++) {
//                            float z = (float) datasetZScores[p];
//                            if (currentWP.getFlipSNPAlleles()[d]) {
//                                z *= -1;
//                            }
//                            // System.out.println(p + "\t" + alleleassessed + "\t" + m_probeList[p] + "\t" + z + "\t" + currentWP.getFlipSNPAlleles()[d]);
//                            if (m_cisOnly) {
//                                // take into account that not all probes have been tested..
//                                int[] probes = currentWP.getProbes();
//                                int probeId = probes[p];
//                                String probeName = m_probeList[probeId];
//                                outfile.writeFloat(z);
//                                if (sb == null) {
//                                    sb = new StringBuilder();
//                                } else {
//                                    sb.append("\t");
//                                }
//                                sb.append(probeName);
//                            } else {
//                                outfile.writeFloat(z);
//                            }
//                        }
//
//                        StringBuilder buffer = new StringBuilder();
//                        buffer.append(snpname)
//                                .append("\t")
//                                .append(BaseAnnot.getAllelesDescription(alleles))
//                                .append("\t")
//                                .append(BaseAnnot.toString(minorAllele))
//                                .append("\t")
//                                .append(BaseAnnot.toString(alleleassessed))
//                                .append("\t")
//                                .append(datasetSNP.getNrCalled())
//                                .append("\t")
//                                .append(maf)
//                                .append("\t")
//                                .append(hwe)
//                                .append("\t")
//                                .append(cr)
//                                .append("\t")
//                                .append(datasetZScores.length)
//                                .append("\t");
//                        if (sb != null) {
//                            buffer.append(sb.toString());
//                        } else {
//                            buffer.append("-");
//                        }
//
//                        snpfile.writeln(buffer.toString());
//
//                    }
//                }
//            }
//        }
//    }
}
