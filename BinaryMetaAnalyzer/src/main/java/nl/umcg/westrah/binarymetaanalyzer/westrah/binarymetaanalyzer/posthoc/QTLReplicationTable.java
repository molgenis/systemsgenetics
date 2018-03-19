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
			int nroverlapsignificant = 0;
			for (int f = 0; f < files.length; f++) {
				EQTL other = output[e][f];
				if (other == null) {
					ln += "\t-"
							+ "\t-"
							+ "\t-"
							+ "\t-"
							+ "\t-";
				} else {
					
					if (!(other.getFDR() < fdrthreshold || includenonSignificanteffects)) {
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
								nroverlapsignificant++;
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
				ln += "\t" + nroverlap + "\t" + nroverlapsignificant;
				
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
