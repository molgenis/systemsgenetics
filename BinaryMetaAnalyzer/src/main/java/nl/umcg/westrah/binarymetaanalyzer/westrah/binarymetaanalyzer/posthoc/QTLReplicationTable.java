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
					boolean onlyoutputQTLThatOverlapAllOtherFiles,
					double fdrthreshold) throws IOException {
		
		
		QTLTextFile t = new QTLTextFile(referenceFile, QTLTextFile.R);
		EQTL[] referenceQTLArr = t.read();
		t.close();
		
		// index
		HashMap<String, Integer> snpGeneToQTL = new HashMap<String, Integer>();
		int ctr = 0;
		ArrayList<EQTL> referenceQTL = new ArrayList<EQTL>();
		for (EQTL e : referenceQTLArr) {
			if (e.getFDR() < fdrthreshold) {
				snpGeneToQTL.put(e.getRsName() + "_" + e.getProbe(), ctr);
				referenceQTL.add(e);
				ctr++;
			}
		}
		
		DetermineLD d = new DetermineLD();
		
		String[] files = otherFiles.split(";");
		String[] filenames = otherFilesNames.split(";");
		EQTL[][] output = new EQTL[referenceQTL.size()][files.length];
		
		for (int f = 0; f < files.length; f++) {
			QTLTextFile t2 = new QTLTextFile(files[f], QTLTextFile.R);
			EQTL[] qtl2 = t2.read();
			t2.close();
			
			for (EQTL e : qtl2) {
				if (e.getFDR() < fdrthreshold) {
					String query = e.getRsName() + "_" + e.getProbe();
					Integer id = snpGeneToQTL.get(query);
					if (id != null) {
						output[id][f] = e;
					}
				}
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
					+ "\t";
			header += "\tZ" +
					"\tRSq" +
					"\tP" +
					"\tFDR";
		}
		
		TextFile out = new TextFile(outputloc, TextFile.W);
		out.writeln(header2);
		out.writeln(header);
		
		for (int e = 0; e < output.length; e++) {
			
			EQTL reference = referenceQTL.get(e);
			
			int n = Descriptives.sum(Primitives.toPrimitiveArr(reference.getDatasetsSamples()));
			
			String ln = reference.getProbe()
					+ "\t" + reference.getProbeChr()
					+ "\t" + reference.getProbeChrPos()
					+ "\t" + reference.getRsName()
					+ "\t" + reference.getRsChr()
					+ "\t" + reference.getRsChrPos()
					+ "\t" + reference.getAlleles()
					+ "\t" + reference.getAlleleAssessed()
					+ "\t" + reference.getZscore()
					+ "\t" + ZScores.zToR(reference.getZscore(), n)
					+ "\t" + reference.getPvalue()
					+ "\t" + reference.getFDR();
			
			for (int f = 0; f < files.length; f++) {
				EQTL other = output[e][f];
				if (other == null) {
					ln += "\t-"
							+ "\t-"
							+ "\t-"
							+ "\t-";
				} else {
					Boolean flip = BaseAnnot.flipalleles(reference.getAlleles(), reference.getAlleleAssessed(), other.getAlleles(), other.getAlleleAssessed());
					if (flip != null) {
						int nother = Descriptives.sum(Primitives.toPrimitiveArr(reference.getDatasetsSamples()));
						double z = other.getZscore();
						if (flip) {
							z *= -1;
						}
						ln += "\t" + z
								+ "\t" + ZScores.zToR(z, nother)
								+ "\t" + other.getPvalue()
								+ "\t" + other.getFDR();
					}
				}
			}
			out.writeln(ln);
		}
		out.close();
	}
	
}
