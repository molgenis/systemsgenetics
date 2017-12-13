package eqtlmappingpipeline.util;


import eqtlmappingpipeline.binarymeta.meta.graphics.ZScorePlot;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;

import java.io.IOException;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class EQTMGeneLinker {
	
	public static void main(String[] args) {
		
		
		EQTMGeneLinker l = new EQTMGeneLinker();
		String eqtl = "C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\eqtl\\eQTLsFDR0.05-PrunedLevel_2CohortFilter_ParalogueFilter_20171204.txt.gz";
		String eqtm = "C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\eQTM\\cis-eQTM250k\\eQTLsFDR0.05-ProbeLevel.txt";
//		String outlink = "C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\eQTM\\linksWithEQTL250k\\";
		String outlink = "C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\eQTM\\linksWithEQTL\\";
		
		String meqtl = "C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\trans250k\\eQTLsFDR-ProbeLevel.txt.gz";
		String compout = "C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\trans250k\\eqtlcompare\\All";
		
		try {
//			l.link(eqtl, eqtm, outlink);
			
			
			l.compare2(eqtl, eqtm, meqtl, compout);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private void compare2(String eqtlfile, String eqtmfile, String meqtlfile, String compout) throws Exception {
		
		
		QTLTextFile f = new QTLTextFile(eqtlfile, QTLTextFile.R);
		ArrayList<EQTL> eqtl = f.readList();
		f.close();
		System.out.println(eqtl.size() + " eQTL ");
		
		f = new QTLTextFile(eqtmfile, QTLTextFile.R);
		ArrayList<EQTL> eqtm = f.readList();
		f.close();
		System.out.println(eqtm.size() + " eQTM ");
		
		f = new QTLTextFile(meqtlfile, QTLTextFile.R);
		ArrayList<EQTL> meqtl = f.readList();
		f.close();
		
		System.out.println(meqtl.size() + " meQTL");
		
		// link meQTL toGenes
		ZScorePlot zp = new ZScorePlot();
		String plotname = compout + "plot.jpg";
		zp.init(3, new String[]{"eQTL", "eQTM", "meQTL"}, false, plotname);
		
		TextFile out = new TextFile(compout + "comparisons.txt", TextFile.W);
		
		String header = "snp" +
				"\tsnpChr" +
				"\tsnpPos" +
				"\tgene" +
				"\tgeneChr" +
				"\tgenePos" +
				"\tcg" +
				"\tcgChr" +
				"\tcgPos" +
				"\tAlleles,assessed(eQTL)" +
				"\teqtlZ" +
				"\teqtmZ" +
				"\tAlleles,assessed(meQTL)" +
				"\tmeqtlZ" +
				"\tmeqtlZFlipped(basedOnAlleles)" +
				"\tmeqtlZFlipped(basedOnAllelesAndEqtm)" +
				"\teqtlFDR" +
				"\teqtmFDR" +
				"\tmeqtlFDR";
		out.writeln(header);
		
		int q = 0;
		int written = 0;
		for (EQTL e : eqtl) {
			
			// link to eQTM
			EQTL maxeqtm = null;
			double minp = 2;
			for (EQTL m : eqtm) {
				if (e.getProbe().equals(m.getProbe())) {
					if (m.getPvalue() < minp) {
						minp = m.getPvalue();
						maxeqtm = m;
					}
				}
			}
			
			EQTL maxmeqtl = null;
			if (maxeqtm != null) {
				double minmeqtlpval = 2;
				for (EQTL m : meqtl) {
					String esnp = e.getRsName();
					String mesnp = m.getRsName();
					String meCG = m.getProbe();
					String eqtmCG = maxeqtm.getRsName();
					if (esnp.equals(mesnp)) {
						if (meCG.equals(eqtmCG)) {
							if (m.getPvalue() < minmeqtlpval) {
								minmeqtlpval = m.getPvalue();
								maxmeqtl = m;
							}
						}
					}
				}
			}
			
			
			// prepare output
			if (maxmeqtl != null) {
				// flip z-score depending on alleles
				Boolean flip = BaseAnnot.flipalleles(e.getAlleles(), e.getAlleleAssessed(), maxmeqtl.getAlleles(), maxmeqtl.getAlleleAssessed());
				double zscoreflipped1 = maxmeqtl.getZscore();
				
				if (flip) {
					zscoreflipped1 *= -1;
				}
				double zscoreflipped2 = zscoreflipped1;
				
				// flip based on eqtm
				if (maxeqtm.getZscore() < 0) {
					zscoreflipped2 *= -1;
				}
				
				zp.draw(e.getZscore(), maxeqtm.getZscore(), 0, 1);
				zp.draw(e.getZscore(), zscoreflipped2, 0, 2);
				zp.draw(maxeqtm.getZscore(), zscoreflipped2, 1, 2);
				
				String lnout = e.getRsName()
						+ "\t" + e.getRsChr()
						+ "\t" + e.getRsChrPos()
						+ "\t" + e.getProbe()
						+ "\t" + e.getProbeChr()
						+ "\t" + e.getProbeChrPos()
						+ "\t" + maxmeqtl.getProbe()
						+ "\t" + maxmeqtl.getProbeChr()
						+ "\t" + maxmeqtl.getProbeChrPos()
						+ "\t" + e.getAlleles() + ", " + e.getAlleleAssessed()
						+ "\t" + e.getZscore()
						+ "\t" + maxeqtm.getZscore()
						+ "\t" + maxmeqtl.getAlleles() + ", " + maxmeqtl.getAlleleAssessed()
						+ "\t" + maxmeqtl.getZscore()
						+ "\t" + zscoreflipped1
						+ "\t" + zscoreflipped2
						+ "\t" + e.getFDR()
						+ "\t" + maxeqtm.getFDR()
						+ "\t" + maxmeqtl.getFDR();
				out.writeln(lnout);
				written++;
			}
			q++;
			System.out.print("\r" + q + " out of " + eqtl.size() + "\t" + written + " written");
		}
		System.out.println();
		System.out.println("Done");
		zp.write(plotname);
		
		out.close();
	}
	
	
	// create SNP/cg/gene combinations using significant eQTL and eQTM
	public void link(String eqtlfile, String eqtmfile, String outfile) throws IOException {
		
		// eqtm and eqtl files share same probes/genes, so we'd need to match on that
		
		
		HashMap<String, ArrayList<QTLObj>> geneToQTL = new HashMap<>();
		TextFile tf = new TextFile(eqtlfile, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		int ctr = 0;
		
		while (elems != null) {
			String snp = elems[1];
			String gene = elems[4];
			
			ArrayList<QTLObj> qtl = geneToQTL.get(gene);
			if (qtl == null) {
				qtl = new ArrayList<>();
			}
			
			double p = Double.parseDouble(elems[0]);
			double z = Double.parseDouble(elems[10]);
			double f = 0;
			try {
				f = Double.parseDouble(elems[elems.length - 1]);
			} catch (NumberFormatException e) {
			
			}
			String alleles = elems[8];
			String alleleAssessed = elems[9];
			QTLObj obj = new QTLObj(snp, alleles, alleleAssessed, p, z, f);
			
			qtl.add(obj);
			geneToQTL.put(gene, qtl);
			
			elems = tf.readLineElems(TextFile.tab);
			ctr++;
		}
		tf.close();
		
		TextFile outAll = new TextFile(outfile + "AllCombos.txt", TextFile.W);
		TextFile outTop = new TextFile(outfile + "TopCombos.txt", TextFile.W);
		String header = "SNP\tCG\tGene\tAllelesEQTL\tAlleleAssessedQTL\tPeQTL\tZScoreEQTL\tFDREQTL\tPeQTM\tZScoreEQTM\tFDREQTM";
		outAll.writeln(header);
		outTop.writeln(header);
		HashSet<String> visitedGenes = new HashSet<String>();
		TextFile tf2 = new TextFile(eqtmfile, TextFile.R);
		tf2.readLine();
		elems = tf2.readLineElems(TextFile.tab);
		while (elems != null) {
			
			String cg = elems[1];
			String gene = elems[4];
			
			ArrayList<QTLObj> objs = geneToQTL.get(gene);
			if (objs != null) {
				
				double p = Double.parseDouble(elems[0]);
				double z = Double.parseDouble(elems[10]);
				double f = Double.parseDouble(elems[elems.length - 1]);
//				String alleles = elems[8];
//				String alleleAssessed = elems[9];
//				String alleles = elems[0];
//				String alleleAssessed = elems[0];
				
				for (QTLObj obj : objs) {
					
					String ln = obj.snp
							+ "\t" + cg
							+ "\t" + gene
							+ "\t" + obj.alleles
							+ "\t" + obj.allelesAssessed
							+ "\t" + obj.p
							+ "\t" + obj.z
							+ "\t" + obj.f
							+ "\t" + p
							+ "\t" + z
							+ "\t" + f;
					
					if (!visitedGenes.contains(gene)) {
						outTop.writeln(ln);
					}
					outAll.writeln(ln);
				}
				
				
				visitedGenes.add(gene);
				
				
			}
			elems = tf2.readLineElems(TextFile.tab);
		}
		outAll.close();
		outTop.close();
		
	}
	
	private class QTLObj {
		
		double p;
		double z;
		double f;
		String snp;
		String alleles;
		String allelesAssessed;
		public String cg;
		public String gene;
		
		public QTLObj(String snp, String alleles, String allelesAssessed, double p, double z, double f) {
			this.p = p;
			this.z = z;
			this.f = f;
			this.snp = snp;
			this.alleles = alleles;
			this.allelesAssessed = allelesAssessed;
		}
		
	}
	
}
