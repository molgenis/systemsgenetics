/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.binarymeta.meta;

import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import eqtlmappingpipeline.binarymeta.meta.graphics.ZScorePlot;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;
import umcg.genetica.io.trityper.bin.BinaryResultSNP;
import umcg.genetica.math.stats.Descriptives;

import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.bin.BinaryResultDataset;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.stats.FisherExactTest;

/**
 *
 * @author harmjan
 */
public class MetaAnalysisPlotThread extends MetaAnalysisCalculationThread {

    public MetaAnalysisPlotThread(LinkedBlockingQueue<MetaAnalysisWorkPackage> input, LinkedBlockingQueue<MetaAnalysisWorkPackage> output,
	    ArrayList<String> snps, ArrayList<String> probes,
	    ArrayList<Byte> snpChr, ArrayList<Integer> snpChrPos,
	    BinaryResultDataset[] ds,
	    Integer[][] snpTranslation,
	    Integer[][] probeTranslationLookupTable, ProbeTranslation probeTranslation,
	    MetaSettings m_settings,
	    ZScorePlot zs, PValueThreshold p) {
	super(input, output, snps, probes, snpChr, snpChrPos, ds, snpTranslation, probeTranslationLookupTable, probeTranslation, m_settings, zs, p);
    }

    @Override
    protected void analyze(MetaAnalysisWorkPackage pack) {

	int s = pack.getSNPNum();

	int[] totalNrSamples = new int[probes.size()];
	double[] zSum = new double[probes.size()];
	double[] zSumAbsolute = new double[probes.size()];
	int[] dsPassQC = new int[probes.size()];
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


	int numDSPassingQC = 0;

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
	}

	boolean[] dsPassQCForProbe = new boolean[ds.length];
	double[][] zScoreForProbe = new double[ds.length][probes.size()];
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
			dsPassQCForProbe[d] = true;
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

					double zscore = zscores[probeId];

					dsPassQC[p]++;

					if (firstSNPPassingQC == null) {
					    firstSNPPassingQC = snpObject;
					} else {
					    Boolean flipalleles = flipalleles(firstSNPPassingQC, snpObject);
					    if (flipalleles == null) {
						System.err.println("ERROR! SNP alleles cannot be matched for snp\t" + snpObject.getName() + "\tin dataset\t" + d);
						System.err.println("This SNP will be excluded from further research");
					    } else if (flipalleles) {
						zscore = -zscore;
						zscoreflipped[d] = true;
					    }
					}

					zScoreForProbe[d][p] = zscore;
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

	pack.clearByteData();

	if (numDSPassingQC > 0) {
	    Double[] metaZPerProbe = null;
	    if (m_settings.isMakezscoretable()) {
		metaZPerProbe = new Double[probes.size()];
	    }
	    for (int p = 0; p < probes.size(); p++) {

		if (dsPassQC[p] > 0 && totalNrSamples[p] > 0) {
		    double zSumVal = zSum[p];
		    double sqrtSample = Descriptives.getSqrt(totalNrSamples[p]);
		    double zScore = zSumVal / sqrtSample;
		    double pValueOverall = Descriptives.convertZscoreToPvalue(zScore);

		    double zSumValAbsolute = zSumAbsolute[p];
		    double zScoreAbs = zSumValAbsolute / sqrtSample;
		    double pValueOverallAbs = Descriptives.convertZscoreToPvalue(zScoreAbs);

		    result[p] = new EQTL();

		    if (pValueOverallAbs < 1) {
			for (int d1 = 0; d1 < ds.length; d1++) {
			    boolean ds1PassesQC = dsPassQCForProbe[d1];
			    if (ds1PassesQC) {
				double zscore = zScoreForProbe[p][d1];
//				if (zscoreflipped[d1]) {
//				    zscore = -zscore;
//				}
				for (int d2 = d1 + 1; d2 < ds.length; d2++) {
				    if (dsPassQCForProbe[d2]) {
					double zscore2 = zScoreForProbe[p][d2];
//					if (zscoreflipped[d2]) {
//					    zscore2 = -zscore2;
//					}
					if (zs != null) {
					    if (pValueOverall < 1E-20) {
						zs.draw(zscore, zscore2, d1, d2);
					    }
					}
				    }
				}
				if (zs != null && pValueOverall < 1E-20) {
				    zs.draw(zscore, zScore, d1, ds.length);
				}
			    }
			}
		    }

		}
	    }

	}

    }
}
