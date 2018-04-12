package nl.umcg.westrah.binarymetaanalyzer.westrah.binarymetaanalyzer.posthoc;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.io.trityper.util.DetermineLD;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class QTLReplicationTable {
	
	
	public void run(String referenceFile,
					String otherFiles,
					String otherFilesNames,
					String outputloc,
					boolean includenonSignificanteffects,
					int minNrDatasetsOverlap,
					double fdrthreshold) throws IOException {
		
		
		System.out.println("Reference: " + referenceFile);
		QTLTextFile t = new QTLTextFile(referenceFile, QTLTextFile.R);
		EQTL[] referenceQTLArr = t.read();
		t.close();
		
		// index
		HashMap<String, Integer> snpGeneToQTL = new HashMap<String, Integer>();
		int ctr = 0;
		ArrayList<EQTL> referenceQTL = new ArrayList<EQTL>();
		for (EQTL e : referenceQTLArr) {
			if (e.getFDR() < fdrthreshold || includenonSignificanteffects) {
				snpGeneToQTL.put(e.getRsName() + "_" + e.getProbe(), ctr);
				referenceQTL.add(e);
				ctr++;
			}
		}
		
		DetermineLD d = new DetermineLD();
		
		String[] files = otherFiles.split(",");
		String[] filenames = otherFilesNames.split(",");
		EQTL[][] output = new EQTL[referenceQTL.size()][files.length];
		
		for (int f = 0; f < files.length; f++) {
			QTLTextFile t2 = new QTLTextFile(files[f], QTLTextFile.R);
			EQTL[] qtl2 = t2.read();
			t2.close();
			
			for (EQTL e : qtl2) {
//				if (e.getFDR() < fdrthreshold || includenonSignificanteffects) {
				String query = e.getRsName() + "_" + e.getProbe();
				Integer id = snpGeneToQTL.get(query);
				if (id != null) {
					output[id][f] = e;
				}
//				}
			}
		}
		
		
		
		/*
		
		rsid chr pos
		gene chr pos
		alleles assessed
		zscore rsq p fdr
		
		per other ds
		zscore rsq p fdr
		
		*/
		
		String header = "Gene" +
				"\tGene-Chr" +
				"\tGene-Pos" +
				"\tRsId" +
				"\tSNP-Chr" +
				"\tSNP-Pos" +
				"\tAlleles" +
				"\tAlleleAssessed" +
				"\tZ" +
				"\tRSq" +
				"\tP" +
				"\tFDR";
		String header2 = "" +
				"\t" +
				"\t" +
				"\t" +
				"\t" +
				"\t" +
				"\t" +
				"\t" +
				"\teQTLGen" +
				"\t" +
				"\t" +
				"\t";
		for (int f = 0; f < files.length; f++) {
			header2 += "\t" + filenames[f]
					+ "\t"
					+ "\t"
					+ "\t"
					+ "\t";
			header += "\tZ" +
					"\tRSq" +
					"\tP(differentEffectSize)" +
					"\tP" +
					"\tFDR";
		}
		header += "\tNrDatasetsTested";
		header += "\tNrDatasetsWithP(differentEffectSize)<0.05";
		header += "\tNrDatasetsWithFDR<" + fdrthreshold;
		header += "\tNrDatasetsWithFDR<" + fdrthreshold + "AndSameDirection";
		header += "\tNrDatasetsWithSameDirection";
		
		
		TextFile out = new TextFile(outputloc, TextFile.W);
		out.writeln(header2);
		out.writeln(header);
		
		ArrayList<ArrayList<Double>> allzscores = new ArrayList<>();
		for (int f = 0; f < filenames.length; f++) {
			allzscores.add(new ArrayList<>());
		}
		
		for (int e = 0; e < output.length; e++) {
			
			EQTL reference = referenceQTL.get(e);
			
			int n = Descriptives.sum(Primitives.toPrimitiveArr(reference.getDatasetsSamples()));
			double r = ZScores.zToR(reference.getZscore(), n);
			String ln = reference.getProbe()
					+ "\t" + reference.getProbeChr()
					+ "\t" + reference.getProbeChrPos()
					+ "\t" + reference.getRsName()
					+ "\t" + reference.getRsChr()
					+ "\t" + reference.getRsChrPos()
					+ "\t" + reference.getAlleles()
					+ "\t" + reference.getAlleleAssessed()
					+ "\t" + reference.getZscore()
					+ "\t" + (r * r)
					+ "\t" + reference.getPvalue()
					+ "\t" + reference.getFDR();
			
			int nroverlap = 0;
			int nrSignificantDifferentEffectSize = 0;
			int nrSignificantFDR = 0;
			int nrSignificantFDRSameDirection = 0;
			int nrSameDirection = 0;
			for (int f = 0; f < files.length; f++) {
				EQTL other = output[e][f];
				if (other == null) {
					ln += "\t-"
							+ "\t-"
							+ "\t-"
							+ "\t-"
							+ "\t-";
				} else {
					
					
					if (!(other.getFDR() >= fdrthreshold || includenonSignificanteffects)) {
						ln += "\tNOTSIGNIFICANT"
								+ "\tNOTSIGNIFICANT"
								+ "\tNOTSIGNIFICANT"
								+ "\tNOTSIGNIFICANT"
								+ "\tNOTSIGNIFICANT";
					} else {
						
						
						nroverlap++;
						Boolean flip = BaseAnnot.flipalleles(reference.getAlleles(), reference.getAlleleAssessed(), other.getAlleles(), other.getAlleleAssessed());
						if (flip != null) {
							
							
							int nother = Descriptives.sum(Primitives.toPrimitiveArr(reference.getDatasetsSamples()));
							double z = other.getZscore();
							if (flip) {
								z *= -1;
							}
							
							
							ArrayList<Double> zscoresforeqtlfile = allzscores.get(f);
							zscoresforeqtlfile.add(z * z);
							double rother = ZScores.zToR(z, nother);
							
							double rzref = zFromCorr(r);
							double rzother = zFromCorr(rother);
							double se = Math.sqrt((1 / (n - 3d)) + (1 / (nother - 3d)));
							double zDiff = (rzother - rzref) / se;
							double rp = cern.jet.stat.Probability.normal(-Math.abs(zDiff)) * 2d;
							
							if (rp < 0.05) {
								nrSignificantDifferentEffectSize++;
							}
							
							ln += "\t" + z
									+ "\t" + (rother * rother)
									+ "\t" + rp
									+ "\t" + other.getPvalue()
									+ "\t" + other.getFDR();
							boolean samedirection = false;
							if ((z > 0 && reference.getZscore() > 0) || (z < 0 && reference.getZscore() < 0)) {
								samedirection = true;
								nrSameDirection++;
							}
							if (other.getFDR() < fdrthreshold) {
								nrSignificantFDR++;
								if (samedirection) {
									nrSignificantFDRSameDirection++;
								}
							}
						} else {
							ln += "\tINCOMPATIBLEALLELES"
									+ "\tINCOMPATIBLEALLELES"
									+ "\tINCOMPATIBLEALLELES"
									+ "\tINCOMPATIBLEALLELES"
									+ "\tINCOMPATIBLEALLELES";
						}
					}
				}
			}
			if (nroverlap >= minNrDatasetsOverlap) {
				ln += "\t" + nroverlap + "\t" + nrSignificantDifferentEffectSize + "\t" + nrSignificantFDR + "\t" + nrSignificantFDRSameDirection + "\t" + nrSameDirection;
				
				out.writeln(ln);
			}
		}
		out.close();
		
		
		// determine lambdas
		TextFile outf2 = new TextFile(outputloc + "-lambdas.txt", TextFile.W);
		outf2.writeln("Tissue/CellType\tNrEQTLs\tLambda(MedianChiSquared)");
		for (int f = 0; f < filenames.length; f++) {
			ArrayList<Double> zscoresforeqtlfile = allzscores.get(f);
			
			double[] z = Primitives.toPrimitiveArr(zscoresforeqtlfile);
			
			double zmed = JSci.maths.ArrayMath.median(z);
			String ln = filenames[f] + "\t" + zscoresforeqtlfile.size() + "\t" + zmed;
			outf2.writeln(ln);
		}
		
		outf2.close();
		
	}
	
	public void rungtex(String referenceFile,
						String otherFiles,
						String otherFilesNames,
						String outputloc,
						boolean includenonSignificanteffects,
						int minNrDatasetsOverlap,
						double pvaluethreshold) throws IOException {
		
		System.out.println("Reference: " + referenceFile);
		QTLTextFile t = new QTLTextFile(referenceFile, QTLTextFile.R);
		EQTL[] referenceQTLArr = t.read();
		t.close();
		System.out.println(referenceQTLArr.length + " eQTLs... ");
		
		// index
		HashMap<String, Integer> snpGeneToQTL = new HashMap<String, Integer>();
		int ctr = 0;
		ArrayList<EQTL> referenceQTL = new ArrayList<EQTL>();
		for (EQTL e : referenceQTLArr) {
			snpGeneToQTL.put(e.getRsName() + "_" + e.getProbe(), ctr);
			referenceQTL.add(e);
			ctr++;
		}
		
		DetermineLD d = new DetermineLD();
		
		String[] files = otherFiles.split(",");
		String[] filenames = otherFilesNames.split(",");
		EQTL[][] output = new EQTL[referenceQTL.size()][files.length];
		
		for (int f = 0; f < files.length; f++) {
			
			QTLTextFile t2 = new QTLTextFile(files[f], QTLTextFile.R);
			EQTL[] qtl2 = t2.read();
			t2.close();
			System.out.println(filenames[f] + "\t" + files[f] + "\t" + qtl2.length + " eQTLs");
			for (EQTL e : qtl2) {
//				if (e.getFDR() < fdrthreshold || includenonSignificanteffects) {
				String query = e.getRsName() + "_" + e.getProbe();
				Integer id = snpGeneToQTL.get(query);
				if (id != null) {
					output[id][f] = e;
				}
//				}
			}
		}
		
		
		
		
		/*
		
		rsid chr pos
		gene chr pos
		alleles assessed
		zscore rsq p fdr
		
		per other ds
		zscore rsq p fdr
		
		*/
		
		String header = "Gene" +
				"\tGene-Chr" +
				"\tGene-Pos" +
				"\tRsId" +
				"\tSNP-Chr" +
				"\tSNP-Pos" +
				"\tAlleles" +
				"\tAlleleAssessed" +
				"\tZ" +
				"\tRSq" +
				"\tP" +
				"\tFDR";
		String header2 = "" +
				"\t" +
				"\t" +
				"\t" +
				"\t" +
				"\t" +
				"\t" +
				"\t" +
				"\teQTLGen" +
				"\t" +
				"\t" +
				"\t";
		for (int f = 0; f < files.length; f++) {
			header2 += "\t" + filenames[f]
					+ "\t"
					+ "\t"
					+ "\t"
					+ "\t";
			header += "\tZ" +
					"\tRSq" +
					"\tP(differentEffectSize)" +
					"\tP" +
					"\tFDR";
		}
		header += "\tNrDatasetsTested";
		header += "\tNrDatsetsSignificantP";
		header += "\tNrDatasetsWithP(differentEffectSize)<0.05";
		
		TextFile out = new TextFile(outputloc, TextFile.W);
		out.writeln(header2);
		out.writeln(header);
		
		for (int e = 0; e < output.length; e++) {
			
			EQTL reference = referenceQTL.get(e);
			
			int n = Descriptives.sum(Primitives.toPrimitiveArr(reference.getDatasetsSamples()));
			double r = ZScores.zToR(reference.getZscore(), n);
			String ln = reference.getProbe()
					+ "\t" + reference.getProbeChr()
					+ "\t" + reference.getProbeChrPos()
					+ "\t" + reference.getRsName()
					+ "\t" + reference.getRsChr()
					+ "\t" + reference.getRsChrPos()
					+ "\t" + reference.getAlleles()
					+ "\t" + reference.getAlleleAssessed()
					+ "\t" + reference.getZscore()
					+ "\t" + (r * r)
					+ "\t" + reference.getPvalue()
					+ "\t" + reference.getFDR();
			
			int nroverlap = 0;
			int nroverlapsignificantlyDifferentEffectSize = 0;
			int nrOverlapSignificant = 0;
			for (int f = 0; f < files.length; f++) {
				EQTL other = output[e][f];
				if (other == null) {
					ln += "\t-"
							+ "\t-"
							+ "\t-"
							+ "\t-"
							+ "\t-";
				} else {
					
					
					if (!(other.getPvalue() < pvaluethreshold || includenonSignificanteffects)) {
						// if pvalue > threshold
						// or don't include significant effects
						ln += "\tNS"
								+ "\tNS"
								+ "\tNS"
								+ "\tNS"
								+ "\tNS";
					} else {
						nroverlap++;
						Boolean flip = BaseAnnot.flipalleles(reference.getAlleles(), reference.getAlleleAssessed(), other.getAlleles(), other.getAlleleAssessed());
						if (flip != null) {
							int nother = Descriptives.sum(Primitives.toPrimitiveArr(reference.getDatasetsSamples()));
							double z = other.getZscore();
							if (flip) {
								z *= -1;
							}
							
							double rother = ZScores.zToR(z, nother);
							
							double rzref = zFromCorr(r);
							double rzother = zFromCorr(rother);
							double se = Math.sqrt((1 / (n - 3d)) + (1 / (nother - 3d)));
							double zDiff = (rzother - rzref) / se;
							double rp = cern.jet.stat.Probability.normal(-Math.abs(zDiff)) * 2d;
							
							if (rp < 0.05) {
								nroverlapsignificantlyDifferentEffectSize++;
							}
							
							if (other.getPvalue() < pvaluethreshold) {
								nrOverlapSignificant++;
							}
							
							ln += "\t" + z
									+ "\t" + rother
									+ "\t" + (rp * rp)
									+ "\t" + other.getPvalue()
									+ "\t" + other.getFDR();
						}
					}
				}
			}
			if (nroverlap >= minNrDatasetsOverlap) {
				ln += "\t" + nroverlap + "\t" + nrOverlapSignificant + "\t" + nroverlapsignificantlyDifferentEffectSize;
				
				out.writeln(ln);
			}
		}
		out.close();
	}
	
	
	private double zFromCorr(double corr1) {
		double raplus = 1 * corr1 + 1;
		double raminus = 1 - corr1;
		double z = (Math.log(raplus) - Math.log(raminus)) / 2;
		return z;
	}
}
