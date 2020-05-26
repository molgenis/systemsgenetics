/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.binarymeta.meta;

import eqtlmappingpipeline.binarymeta.meta.graphics.ZScorePlot;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.bin.BinaryResultDataset;
import umcg.genetica.io.trityper.bin.BinaryResultSNP;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.stats.Descriptives;

import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;

/**
 * @author harmjan
 */
public class MetaAnalysisCalculationThread extends Thread {

	protected LinkedBlockingQueue<MetaAnalysisWorkPackage> m_queue_input;
	protected LinkedBlockingQueue<MetaAnalysisWorkPackage> m_queue_output;
	protected ArrayList<String> probes;
	protected ArrayList<String> snps;
	protected ArrayList<Byte> snpChr;
	protected ArrayList<Integer> snpChrPos;
	protected BinaryResultDataset[] ds;
	protected Integer[][] snpTranslation;
	protected Integer[][] probeTranslationLookupTable;
	protected ProbeTranslation probeTranslation;
	protected MetaSettings m_settings;
	protected ZScorePlot zs;
	protected Inflater inflater = new Inflater();
	protected PValueThreshold pvaluethreshold;
	private int numEffects = 0;
	private int numSNPs = 0;

	public MetaAnalysisCalculationThread(LinkedBlockingQueue<MetaAnalysisWorkPackage> input, LinkedBlockingQueue<MetaAnalysisWorkPackage> output,
										 ArrayList<String> snps, ArrayList<String> probes,
										 ArrayList<Byte> snpChr, ArrayList<Integer> snpChrPos,
										 BinaryResultDataset[] ds,
										 Integer[][] snpTranslation,
										 Integer[][] probeTranslationLookupTable, ProbeTranslation probeTranslation,
										 MetaSettings m_settings,
										 ZScorePlot zs, PValueThreshold p) {
		this.probes = probes;
		this.snps = snps;
		this.snpChr = snpChr;
		this.snpChrPos = snpChrPos;
		this.snpTranslation = snpTranslation;
		this.ds = ds;
		this.probeTranslation = probeTranslation;
		this.probeTranslationLookupTable = probeTranslationLookupTable;
		this.m_settings = m_settings;
		this.zs = zs;
		this.pvaluethreshold = p;
		m_queue_input = input;
		m_queue_output = output;
	}

	@Override
	public void run() {
		boolean poison = false;
		while (!poison) {
			try {
				MetaAnalysisWorkPackage pack = m_queue_input.take();
				poison = pack.getPoison();
				if (!poison) {
					analyze(pack);

//                    if(taken % printperiterations == 0){
//                        System.out.println("Thread "+this.getName()+" calculated "+taken+" workpackages.");
//                    }
				}

			} catch (InterruptedException ex) {
				ex.printStackTrace();
			}
		}

		System.out.println(this.getName() + " - Poisoned - Num tests passed QC: " + numEffects + "\t" + numSNPs);
	}

	protected void analyze(MetaAnalysisWorkPackage pack) {

		int s = pack.getSNPNum();


		// DEBUG
//	boolean verbose = false;
//	if (snps.get(s).equals("rs6919346")) {
//	    verbose = true;
//	}


		int[] totalNrSamples = new int[probes.size()];
		double[] zSum = new double[probes.size()];
		double[] zSumAbsolute = new double[probes.size()];
		int[] dsPassQC = new int[probes.size()];
		Result r = new Result();
		r.finalzscores = new Double[probes.size()];
		r.finalpvalues = new Double[probes.size()];
		r.numSamples = new Integer[probes.size()][ds.length];
		r.datasetZScores = new Double[probes.size()][ds.length];
		r.dspassingqc = new boolean[probes.size()][ds.length];
		r.snp = s;
		r.passesQC = true;
		r.datasets = new String[ds.length];
		boolean[] zscoreflipped = new boolean[ds.length];
		EQTL[] result = new EQTL[probes.size()];

		BinaryResultSNP firstSNPPassingQC = null;

		Byte snpchr = snpChr.get(s);
		Integer snpchrpos = snpChrPos.get(s);
		boolean snphaspropermapping = true;
		if (snpchr == null || snpchrpos == null || snpchr == -1) {
			snpchr = -1;
			snpchrpos = -1;
			snphaspropermapping = false;
		}

		StringBuilder zscoretableout = new StringBuilder();

		int numDSPassingQC = 0;

		HashSet<Integer> probesTestedHash = new HashSet<Integer>();
		boolean[] testprobes = new boolean[probes.size()];
		for (int p = 0; p < probes.size(); p++) {
			byte probechr = probeTranslation.getProbeChr(p);
			int probechrpos = probeTranslation.getProbeChrPos(p);
			boolean testprobe = false;

			if (m_settings.isCis() && m_settings.isTrans()) {
				testprobe = true;
			} else if (m_settings.isCis() && !m_settings.isTrans()) {
				if (snpchr < 1 || probechr < 1) {
					testprobe = false;
				} else if (probechr == snpchr) {
					if (Math.abs(snpchrpos - probechrpos) < m_settings.getCisdistance()) {
						testprobe = true;
					} else {
						testprobe = false;
					}
				} else {
					testprobe = false;
				}
			} else if (!m_settings.isCis() && m_settings.isTrans()) {
				if (snpchr < 1 || probechr < 1) {
					testprobe = false;
				} else if (probechr == snpchr) {
					if (Math.abs(snpchrpos - probechrpos) > m_settings.getTransdistance()) {
						testprobe = true;
					} else {
						testprobe = false;
					}
				} else {
					testprobe = true;
				}
			}
			testprobes[p] = testprobe;

			if (testprobe) {
				probesTestedHash.add(p);
			}
		}


		for (int d = 0; d < ds.length; d++) {

			Integer snpId = snpTranslation[d][s];

			if (snpId != null) {

				BinaryResultSNP snpObject = pack.getSNPObject(d); // ds[d].getSnps()[snpId];

//		long pointer = snpObject.getzScoreIndex();
//		long nextpointer = -1;
//
//		if (snpId + 1 < ds[d].getSnps().length) {
//		    SNP snpObject2 = ds[d].getSnps()[snpId + 1];
//		    nextpointer = snpObject2.getzScoreIndex();
//		}

				byte[] data = pack.getData(d);
				Float[] zscores = null;
				if (data != null) {
					try {
						zscores = inflate(data, ds[d].getNumProbes()); //
						pack.setData(d, null);
					} catch (DataFormatException ex) {
						Logger.getLogger(MetaAnalysisCalculationThread.class.getName()).log(Level.SEVERE, null, ex);
					}


					if (zscores != null) {
						numDSPassingQC++;
						// weight for dataset d
						int nrSamples = snpObject.getNumsamples();
						double weight = Descriptives.getSqrt(nrSamples);

						for (int p = 0; p < probes.size(); p++) {

							boolean testprobe = testprobes[p];
							if (testprobe) {
								Integer probeId = probeTranslationLookupTable[d][p];

								if (!testprobe && probeId != null) {
									zscores[probeId] = null;
								} else if (probeId != null && testprobe) {
									if (zscores[probeId] != null) {

										totalNrSamples[p] += nrSamples;
										r.dspassingqc[p][d] = true;
										r.numSamples[p][d] = nrSamples;

										double zscore = zscores[probeId];


										r.datasets[d] = ds[d].getM_name().intern();
										dsPassQC[p]++;

										if (firstSNPPassingQC == null) {
											firstSNPPassingQC = snpObject;
										} else {
											Boolean flipalleles = flipalleles(firstSNPPassingQC, snpObject);
											if (flipalleles == null) {
												System.err.println("ERROR! SNP alleles cannot be matched for snp\t" + snpObject.getName() + "\tin dataset\t" + d);
												System.err.println("This SNP will be excluded from further research");
												r.passesQC = false;
											} else if (flipalleles) {
												zscore = -zscore;
												zscoreflipped[d] = true;
											}
										}

										r.datasetZScores[p][d] = new Double(zscore);

//					if (verbose) {
//					    System.out.println(d + "\t" + r.datasetZScores[p][d]);
//					}
										zSumAbsolute[p] += Math.abs(zscore * weight);
										zSum[p] += (zscore * weight);
									} else {
									}
								}
							}
						}
						for (int i = 0; i < zscores.length; i++) {
							zscores[i] = null;
						}
					}
				}


			}
		}

//	if (verbose) {
////	    System.exit(0);
//	}

		pack.clearByteData();

		int numDSThatMinimallyShouldHaveEffect = m_settings.getSnpDatasetPresenceThreshold();
		if (numDSThatMinimallyShouldHaveEffect == 0) {
			numDSThatMinimallyShouldHaveEffect = 1;
		}

		if (numDSPassingQC >= numDSThatMinimallyShouldHaveEffect) {
			pack.setPassedQC(true);
			Double[] metaZPerProbe = null;
			if (m_settings.isMakezscoretable()) {
				metaZPerProbe = new Double[probes.size()];
			}
			int probesTested = 0;
			numSNPs++;
			for (int p = 0; p < probes.size(); p++) {

				if (dsPassQC[p] >= numDSThatMinimallyShouldHaveEffect && totalNrSamples[p] > 0) {
					numEffects++;
					probesTestedHash.add(p);
					probesTested++;
					double zSumVal = zSum[p];
					double sqrtSample = Descriptives.getSqrt(totalNrSamples[p]);
					double metaZScore = zSumVal / sqrtSample;
					double pValueOverall = Descriptives.convertZscoreToPvalue(metaZScore);

					double zSumValAbsolute = zSumAbsolute[p];
					double zScoreAbs = zSumValAbsolute / sqrtSample;
					double pValueOverallAbs = Descriptives.convertZscoreToPvalue(zScoreAbs);


					boolean outputeqtl = false;
					if (m_settings.isMakezscoretable()) {
						outputeqtl = true;
					} else if (pValueOverall <= pvaluethreshold.getPvalue()) {
						outputeqtl = true;
					}

					if (outputeqtl) {
						result[p] = new EQTL();
						EQTL e = result[p];
						e.setRsChr(snpChr.get(s));
						e.setRsChrPos(snpChrPos.get(s));
						e.setProbeChr(probeTranslation.getProbeChr(p));
						e.setProbeChrPos(probeTranslation.getProbeChrPos(p));
						e.setDatasets(r.datasets);
						e.setAlleleAssessed(BaseAnnot.toString(firstSNPPassingQC.getAssessedAllele()).intern());
						byte[] alleles = firstSNPPassingQC.getAlleles();
						String alleleStr = (BaseAnnot.toString(alleles[0]) + "/" + BaseAnnot.toString(alleles[1])).intern();
						e.setAlleles(alleleStr);
						e.setDatasetZScores(r.datasetZScores[p]);
						e.setZscore(metaZScore);
						e.setPvalue(pValueOverall);
						e.setZscoreAbs(zScoreAbs);
						e.setPvalueAbs(pValueOverallAbs);

						if (m_settings.isUseAbsoluteZscore()) {
							e.setUseAbsoluteZScore();
						}

						if (pValueOverallAbs < 1) {
							for (int d1 = 0; d1 < ds.length; d1++) {
								boolean ds1PassesQC = r.dspassingqc[p][d1];
								if (ds1PassesQC) {
									double datasetZScore = r.datasetZScores[p][d1];
									if (zscoreflipped[d1]) {
										datasetZScore = -datasetZScore;
									}
									for (int d2 = d1 + 1; d2 < ds.length; d2++) {
										if (r.dspassingqc[p][d2]) {
											double zscore2 = r.datasetZScores[p][d2];
											if (zscoreflipped[d2]) {
												zscore2 = -zscore2;
											}
											if (zs != null) {
//						if ((datasetZScore < -10 && zscore2 > 10) || (datasetZScore > 10 && zscore2 < -10)) {
//						    System.out.println("");
//						    System.out.println("Opposite effect: ");
//						    System.out.println(probeTranslation.getProbes()[p]
//							    + "\t" + probeTranslation.getProbeChr(p)
//							    + "\t" + probeTranslation.getProbeChrPos(p)
//							    + "\t" + probeTranslation.getProbeSymbol(p));
//						    SNP snpobj1 = pack.getSNPObject(d1);
//						    SNP snpobj2 = pack.getSNPObject(d2);
//						    Integer probe1 = probeTranslationLookupTable[d1][p];
//						    Integer probe2 = probeTranslationLookupTable[d2][p];
//						    System.out.println(d1 + "\t" + datasetZScore
//							    + "\t" + snpobj1.getName()
//							    + "\t" + snpobj1.getChr()
//							    + "\t" + snpobj1.getChrpos()
//							    + "\t" + snpobj1.getHwe()
//							    + "\t" + snpobj1.getMaf()
//							    + "\t" + BaseAnnot.toString(snpobj1.getAlleles()[0]) + "/" + BaseAnnot.toString(snpobj1.getAlleles()[1]) + "-" + BaseAnnot.toString(snpobj1.getMinorAllele()) + "-" + BaseAnnot.toString(snpobj1.getAssessedAllele()) + "\t" + zscoreflipped[d1]);
//
//						    System.out.println(d2 + "\t" + zscore2
//							    + "\t" + snpobj2.getName()
//							    + "\t" + snpobj2.getChr()
//							    + "\t" + snpobj2.getChrpos()
//							    + "\t" + snpobj2.getHwe()
//							    + "\t" + snpobj2.getMaf()
//							    + "\t" + BaseAnnot.toString(snpobj2.getAlleles()[0]) + "/" + BaseAnnot.toString(snpobj2.getAlleles()[1]) + "-" + BaseAnnot.toString(snpobj2.getMinorAllele()) + "-" + BaseAnnot.toString(snpobj2.getAssessedAllele()) + "\t" + zscoreflipped[d2]);
//
//
//						    double maf1 = snpobj1.getMaf();
//						    int nrSamples1 = snpobj1.getNumsamples();
//						    int nrSamplesA1 = (int) Math.ceil(nrSamples1 * maf1);
//						    double maf2 = snpobj2.getMaf();
//						    int nrSamples2 = snpobj2.getNumsamples();
//						    int nrSamplesA2 = (int) Math.ceil(nrSamples2 * maf2);
//
//						    int ctA1 = nrSamplesA1;
//						    int ctB1 = nrSamples1 - nrSamplesA1;
//						    int ctA2 = nrSamplesA2;
//						    int ctB2 = nrSamples2 - nrSamplesA2;
//						    FisherExactTest fet = new FisherExactTest();
//						    double pval = fet.getFisherPValue(ctA1, ctA2, ctB1, ctB2);
//						    System.out.println("Fisher Exact P-value for SNPs:\t" + pval);
//
//						    System.out.println("");
//						}
												if (pValueOverall < 1E-15) {
													zs.draw(new Double(datasetZScore), new Double(zscore2), d1, d2);
												}
											}
										}
									}
									if (zs != null && pValueOverall < 1E-15) {
										zs.draw(new Double(datasetZScore), new Double(metaZScore), d1, ds.length);
									}
								}
							}
						}
						//
						e.setDatasetsSamples(r.numSamples[p]);
						e.setProbe(probes.get(p).intern());
						e.setRsName(firstSNPPassingQC.getName().intern());
						e.setProbeHUGO(probeTranslation.getProbeSymbol(p).intern());

					}

					if (m_settings.isMakezscoretable()) {
						metaZPerProbe[p] = metaZScore;
					}
				} else {
					r.finalzscores[p] = null;
				}
			}

			if (m_settings.isMakezscoretable()) {

				if (firstSNPPassingQC != null) {
					zscoretableout.append(snps.get(s));
					zscoretableout.append("\t").append(BaseAnnot.toString(firstSNPPassingQC.getAlleles()[0])).append("/").append(BaseAnnot.toString(firstSNPPassingQC.getAlleles()[1])).append("\t").append(BaseAnnot.toString(firstSNPPassingQC.getAssessedAllele()));

					for (int i = 0; i < metaZPerProbe.length; i++) {
						zscoretableout.append("\t").append(metaZPerProbe[i]);
						metaZPerProbe[i] = null;
					}
					metaZPerProbe = null;
					pack.setZScoreOut(zscoretableout.toString());
				}
			}
			r.clearData();


			if (numDSPassingQC > 0) {
				pack.setProbesTestedHash(probesTestedHash);
			} else {
				pack.setProbesTestedHash(new HashSet<Integer>());
			}
			pack.setNumOfTestedProbes(probesTested);
			pack.setResult(result);
			try {
				m_queue_output.put(pack);
			} catch (InterruptedException ex) {
				ex.printStackTrace();
			}

		}


	}

	// TODO: AT / GC SNPs??
	public Boolean flipalleles(BinaryResultSNP firstSNPPassingQC, BinaryResultSNP snpObject) {
		byte[] allelesfirst = firstSNPPassingQC.getAlleles();
		byte allelefirstassessed = firstSNPPassingQC.getAssessedAllele();

		byte[] allelessecond = snpObject.getAlleles();
		byte allelesecondassessed = snpObject.getAssessedAllele();

		int nridenticalalleles = 0;

		for (int i = 0; i < allelesfirst.length; i++) {
			byte allele1 = allelesfirst[i];
			for (int j = 0; j < allelessecond.length; j++) {
				if (allelessecond[j] == allele1) {
					nridenticalalleles++;
				}
			}
		}

		if (nridenticalalleles == 2) {
			// alleles are identical. check if same allele was assessed...
			if (allelefirstassessed == allelesecondassessed) {
				return false;
			} else {
				return true;
			}
		} else {
			// try complement
			allelessecond = convertToComplementaryAlleles(allelessecond);
			allelesecondassessed = BaseAnnot.getComplement(allelesecondassessed);
			nridenticalalleles = 0;

			for (int i = 0; i < allelesfirst.length; i++) {
				byte allele1 = allelesfirst[i];
				for (int j = 0; j < allelessecond.length; j++) {
					if (allelessecond[j] == allele1) {
						nridenticalalleles++;
					}
				}
			}

			if (nridenticalalleles == 2) {
				// alleles are identical. check if same allele was assessed...
				if (allelefirstassessed == allelesecondassessed) {
					return false;
				} else {
					return true;
				}
			}
		}
		return null;
	}

	public byte[] convertToComplementaryAlleles(byte[] allelesToCompare) {
		byte[] allelesComplementary = new byte[2];
		for (int a = 0; a < 2; a++) {
			allelesComplementary[a] = BaseAnnot.getComplement(allelesToCompare[a]);
		}
		return allelesComplementary;
	}

	protected Float[] inflate(byte[] buffer, int numElems) throws DataFormatException {
		inflater.setInput(buffer);
		inflater.finished();
		byte[] decompressed = new byte[numElems * 4];
		inflater.inflate(decompressed);

		long actuallydecompressed = inflater.getBytesWritten();
		if (actuallydecompressed != numElems * 4) {
			throw new DataFormatException("IO Error: uncompressed data does not correspond to the size requested\t" + actuallydecompressed + "\t" + numElems * 4);
		}

		inflater.reset();

		ByteBuffer bytebuffer = ByteBuffer.wrap(decompressed);
		Float[] output = new Float[numElems];
		int ctr = 0;
		for (int i = 0; i < numElems; i++) {
			Float f = bytebuffer.getFloat();
			if (f.isNaN()) {
				f = null;
			} else {
				ctr++;
			}
			output[i] = f;
		}

		decompressed = null;

		if (ctr == 0) {
			return null;
		} else {
			return output;
		}
	}
}
