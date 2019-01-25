/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.binarymeta.meta;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.LinkedBlockingQueue;

/**
 * @author harmjan
 */
public class MetaAnalysisResultThread extends Thread {

	private final LinkedBlockingQueue<MetaAnalysisWorkPackage> m_queue_input;
//    private double pvaluethreshold = 1;


	private static String header = "PValue\t"
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

	private int ctr = 0;
	private EQTL[] eQTLBuffer = new EQTL[100000];
	private EQTL[] finalEQTLBuffer = new EQTL[0];
	private int nrInFinalBuffer = 0;
	private static MetaSettings m_settings;
	private int perm;
	private String[] datasets;
	private final TextFile zscoretable;
	private final PValueThreshold pvaluethreshold;
	private final ArrayList<String> snps;
	private final HashMap<String, HashSet<String>> snpProbeSelection;
	private final ArrayList<String> probes;

	public MetaAnalysisResultThread(LinkedBlockingQueue<MetaAnalysisWorkPackage> input,
									MetaSettings m_settings,
									String[] datasets,
									int perm,
									TextFile zscoretable, PValueThreshold p, ArrayList<String> snps, HashMap<String, HashSet<String>> snpProbeSelection, ArrayList<String> probes) {
		this.m_settings = m_settings;
		this.datasets = datasets;
		this.perm = perm;
		this.zscoretable = zscoretable;
		this.pvaluethreshold = p;
		m_queue_input = input;
		this.snps = snps;
		this.snpProbeSelection = snpProbeSelection;
		this.probes = probes;
	}

	TextFile snpout = null;

	@Override
	public void run() {
		boolean poison = false;
		try {
			snpout = new TextFile(m_settings.getOutput() + "snpsandnreqtls.txt", TextFile.W);
			while (!poison) {
				try {
					MetaAnalysisWorkPackage pack = m_queue_input.take();
					if (!pack.getPoison()) {
						Integer snpnum = pack.getSNPNum();
						String snp = snps.get(snpnum);
						if (snpProbeSelection == null || snpProbeSelection.containsKey(snp)) {
							analyze(pack);
						}


//                    if(taken % printperiterations == 0){
//                        System.out.println("Thread "+this.getName()+" calculated "+taken+" workpackages.");
//                    }
					} else {
						poison = pack.getPoison();
//                    System.out.println("Thread " + m_name + " got killed by a poisonous workpackage, but was bravely able to perform\t" + testsPerformed + "\ttests");
					}

				} catch (InterruptedException ex) {
					ex.printStackTrace();
				}
			}

			if (ctr > 0) {
				mergebuffers(ctr);
			}
			snpout.close();

			java.util.Arrays.sort(finalEQTLBuffer);

			// write eQTL results..

			writeresults();

			TextFile out = new TextFile(m_settings.getOutput() + "/NumberOfEQTLSTotal.txt", TextFile.W);
			out.writeln("Number of eQTLs in total: " + totalNumberOfEQTLs);
			System.out.println("Number of eQTLs in total: " + totalNumberOfEQTLs);
			out.writeln("Number of snps in total: " + uniqueSNPs.size());
			out.writeln("Number of snps in total not passing QC: " + uniqueSNPsNotPassingQC.size());


			System.out.println("Number of snps in total: " + uniqueSNPs.size());
			TextFile out2 = new TextFile(m_settings.getOutput() + "/TestedSNPs.txt", TextFile.W);
			List<String> list = new ArrayList<String>(uniqueSNPs);
			out2.writeList(list);
			out2.close();

			out2 = new TextFile(m_settings.getOutput() + "/TestedSNPsNPQC.txt", TextFile.W);
			list = new ArrayList<String>(uniqueSNPsNotPassingQC);
			out2.writeList(list);

			out.writeln("Number of probes in total: " + uniqueProbes.size());
			System.out.println("Number of probes in total: " + uniqueProbes.size());
			out2 = new TextFile(m_settings.getOutput() + "/TestedProbes.txt", TextFile.W);
			List<Integer> list2 = new ArrayList<Integer>(uniqueProbes);
			ArrayList<String> list2str = new ArrayList<String>();
			for (Integer i : list2) {
				list2str.add("" + i);
			}

			out2.writeList(list2str);
			out.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private HashSet<String> uniqueSNPs = new HashSet<String>();
	private HashSet<String> uniqueSNPsNotPassingQC = new HashSet<String>();
	private HashSet<Integer> uniqueProbes = new HashSet<Integer>();
	private int totalNumberOfEQTLs = 0;

	private void analyze(MetaAnalysisWorkPackage pack) {

		Integer snpnum = pack.getSNPNum();
		String snp = snps.get(snpnum).intern();

		HashSet<String> allowedProbes = null;
		if (snpProbeSelection != null) {
			allowedProbes = snpProbeSelection.get(snp);
		}

		Integer[] probeList = pack.getListOfTestedProbes();
		for (int i = 0; i < probeList.length; i++) {
			String probe = probes.get(probeList[i]);
			if (probe != null && (allowedProbes == null || allowedProbes.contains(probe))) {
				totalNumberOfEQTLs++;
				uniqueProbes.add(probeList[i]);
			}
		}

		if (pack.getPassedQC()) {
			uniqueSNPs.add(snps.get(snpnum).intern());
		} else {
			uniqueSNPsNotPassingQC.add(snps.get(snpnum).intern());
		}
		if (m_settings.isMakezscoretable() && zscoretable != null) {
			try {
				String zscoreout = pack.getZScoreOut();
				if (zscoreout != null) {
					zscoretable.writeln(zscoreout);
					pack.setZScoreOut(null);
				}
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		EQTL[] finalEQTLs = pack.getResult();

		int nreQTLsForSNP = 0;
		for (int p = 0; p < finalEQTLs.length; p++) {

//	    if (finalEQTLs[p] != null) {
////		uniqueProbes.add(finalEQTLs[p].getProbe());
//	    }
			if (finalEQTLs[p] != null && finalEQTLs[p].getPvalue() <= pvaluethreshold.getPvalue() && (allowedProbes == null || allowedProbes.contains(finalEQTLs[p].getProbe()))) {
				nreQTLsForSNP++;
				// check cis / trans constraints ...
//                        if(finalEQTLs[p].getProbeChr())
//		boolean includeEQTL = true;
//                        if(transAnalysis && !cisAnalysis){
//                            if(finalEQTLs[p].getProbeChr().byteValue() == -1 || finalEQTLs[p].getRsChr().byteValue() == -1){
//                                includeEQTL = false;
//                            }
//                            if(finalEQTLs[p].getProbeChrPos().intValue() == -1 || finalEQTLs[p].getRsChr().byteValue() == -1){
//                                includeEQTL = false;
//                            }
//                            if(finalEQTLs[p].getProbeChr().byteValue() ==  finalEQTLs[p].getRsChr().byteValue()){
//                                if( Math.abs(finalEQTLs[p].getProbeChr().byteValue() - finalEQTLs[p].getRsChr().byteValue()) < cisProbeDistance){
//                                    includeEQTL = false;
//                                }
//                            }
//                        }

//		if (includeEQTL) {
				eQTLBuffer[ctr] = finalEQTLs[p];
				ctr++;
				if (ctr == eQTLBuffer.length) {
					mergebuffers(ctr);
					ctr = 0;
//                                System.out.println("SNPs tested: "+s+"/"+snps.size()+", threshold: "+pvaluethreshold);
				}
//		}
			} else {
				if (finalEQTLs[p] != null) {
					finalEQTLs[p].clearData();
					finalEQTLs[p] = null;
				}
			}
		}
		finalEQTLs = null;
		try {
			snpout.writeln(snps.get(pack.getSNPNum()) + "\t" + nreQTLsForSNP);
		} catch (Exception e) {
			e.printStackTrace();
		}
		pack.clearData();
		pack = null;
	}

	protected void mergebuffers(int ctr) {
		EQTL[] toMerge = null;
		if (ctr < eQTLBuffer.length) {
			toMerge = new EQTL[ctr];
			System.arraycopy(eQTLBuffer, 0, toMerge, 0, ctr);
		} else {
			toMerge = eQTLBuffer;
		}

		EQTL[] tmp = new EQTL[finalEQTLBuffer.length + toMerge.length];
		System.arraycopy(toMerge, 0, tmp, 0, toMerge.length);
		System.arraycopy(finalEQTLBuffer, 0, tmp, toMerge.length, finalEQTLBuffer.length);


		nrInFinalBuffer += toMerge.length;
		if (nrInFinalBuffer < m_settings.getFinalEQTLBufferMaxLength()) {
			finalEQTLBuffer = tmp;
		} else {

			java.util.Arrays.sort(tmp);
			finalEQTLBuffer = new EQTL[m_settings.getFinalEQTLBufferMaxLength()];
//            System.out.println(finalEQTLBuffer.length+"\t"+tmp.length);
			System.arraycopy(tmp, 0, finalEQTLBuffer, 0, m_settings.getFinalEQTLBufferMaxLength());
			nrInFinalBuffer = m_settings.getFinalEQTLBufferMaxLength();
			pvaluethreshold.setPvalue(finalEQTLBuffer[nrInFinalBuffer - 1].getPvalue());

		}
	}

	private void writeresults() throws IOException {


		TextFile out = null;
		if (perm > 0) {
			out = new TextFile(m_settings.getOutput() + "PermutedEQTLsPermutationRound" + perm + ".txt.gz", TextFile.W);
		} else {
			out = new TextFile(m_settings.getOutput() + "eQTLs.txt.gz", TextFile.W);
		}


		out.write(header + "\n");


		for (int i = 0; i < finalEQTLBuffer.length; i++) {
			finalEQTLBuffer[i].setDatasets(datasets);
			out.writeln(finalEQTLBuffer[i].toString());
		}

		out.close();

		TextFile oppositeEffects = null;
		if (perm > 0) {
			oppositeEffects = new TextFile(m_settings.getOutput() + "OppositeEffects-PermutedEQTLsPermutationRound" + perm + ".txt.gz", TextFile.W);
		} else {
			oppositeEffects = new TextFile(m_settings.getOutput() + "OppositeEffects-eQTLs.txt.gz", TextFile.W);
		}

		for (int i = 0; i < finalEQTLBuffer.length; i++) {
			String oppositeEffectIndicator = "";
			double pValueOverall = finalEQTLBuffer[i].getPvalue();
			double pValueAbs = finalEQTLBuffer[i].getPvalueAbs();
			if (pValueAbs < pValueOverall) {
				oppositeEffectIndicator = "OppositeEffect";
				if (pValueAbs <= pValueOverall / 100000) {
					oppositeEffectIndicator = "StrongOppositeEffect";
				}

				oppositeEffects.writeln(oppositeEffectIndicator + "\t" + pValueAbs + "\t" + finalEQTLBuffer[i].getZscoreAbs() + "\t" + finalEQTLBuffer[i].toString());
			}
		}
		oppositeEffects.close();
	}
}
