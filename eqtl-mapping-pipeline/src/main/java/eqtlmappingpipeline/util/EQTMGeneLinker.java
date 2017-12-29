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
		eqtl = "C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\eqtl\\eQTLsFDR0.05-PrunedLevel_2CohortFilter_20170925.txt.gz";
		
		String eqtm = "C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\eQTM\\cis-eQTM250k\\eQTLsFDR0.05-ProbeLevel.txt";
//		String outlink = "C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\eQTM\\linksWithEQTL250k\\";
		String outlink = "C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\eQTM\\linksWithEQTL\\";
		
		String meqtl = "C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\trans250k\\eQTLsFDR0.05-ProbeLevel.txt";
		String compout = "C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\trans250k\\eqtlcompare\\FDR005";
		
		String truepositivesAndNegatives = "C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\tptn\\2017-11-14-TP - TN Set - v2.txt";
		try {
//			l.link(eqtl, eqtm, outlink);
			
			String outdir = "C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\tptn\\";
			meqtl = "C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\tptn\\transtp\\eQTLsFDR0.05-ProbeLevel.txt";
//			l.compareTP(eqtm, meqtl, truepositivesAndNegatives, outdir);
			compout = "C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\tptn\\FDR005";
//			l.compare2(eqtl, eqtm, meqtl, truepositivesAndNegatives, compout);

//			String meprs = "C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\mePRS\\eQTLsFDR0.05.txt";
//			String out = "C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\mePRS\\genes.txt";
//			l.linkmePRSToEQTM(meprs, eqtm, out);
//
//			String ePRS = "C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\mePRS\\ePRS\\eQTLsFDR0.05-ProbeLevel.txt";
			
			String eqtlfile = "C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\eqtl\\lclcompare\\eQTLGenVersusNewLCLMeta-FDR0.05-eQTLsWithIdenticalDirecton.txt.gz";
			String meqtlfile = "C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\trans250k\\eQTLsFDR0.05-ProbeLevel.txt";
			String eqtmfile = "C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\eQTM\\cis-eQTM250k\\eQTLsFDR0.05-ProbeLevel.txt";
			String out = "C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\eqtl\\lclcompare\\meqtl\\2017-12-22-meQTLvsLCLEQTL";
			l.compare2(eqtlfile, eqtmfile, meqtlfile, null, out);
			
			
			meqtlfile = "C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\trans250k\\eQTLsFDR0.05-ProbeLevel.txt";
			eqtmfile = "C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\eQTM\\cis-eQTM250k\\eQTLsFDR0.05-ProbeLevel.txt";
			out = "C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\eqtl\\methylationQTL-FDR0.05.txt";
			
			l.convertToEQTLFile(meqtlfile, eqtmfile, out, 0.05);
			
			
		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private void convertToEQTLFile(String meqtlfile, String eqtmfile, String out, double threshold) throws IOException {
		QTLTextFile f = new QTLTextFile(eqtmfile, QTLTextFile.R);
		
		ArrayList<EQTL> eqtm = f.readList();
		f.close();
		System.out.println(eqtm.size() + " eQTM ");
		
		f = new QTLTextFile(meqtlfile, QTLTextFile.R);
		ArrayList<EQTL> meqtl = f.readList();
		f.close();
		
		
		// match eQTM to meQTL
		TextFile inq = new TextFile(meqtlfile, TextFile.R);
		TextFile outf = new TextFile(out, TextFile.W);
		outf.writeln(inq.readLine());
		inq.close();
		
		HashSet<String> printed = new HashSet<String>();
		for (int q = 0; q < meqtl.size(); q++) {
			EQTL me = meqtl.get(q);
			EQTL matchingeqtm = null;
			for (EQTL me2 : eqtm) {
				if (me2.getRsName().equals(me.getProbe())) {
					// matching CG
					if (matchingeqtm == null) {
						matchingeqtm = me2;
					} else {
						if (me2.getPvalue() < matchingeqtm.getPvalue()) {
							matchingeqtm = me2;
						}
					}
				}
			}
			
			if (matchingeqtm != null) {
				
				me.setProbe(matchingeqtm.getProbe());
				
				String eqtl = me.getRsName() + "_" + me.getProbe();
				if (!printed.contains(eqtl)) {
					// print
					// flip effect depending on eQTM
					double eqtmz = matchingeqtm.getZscore();
					if (!matchingeqtm.getAlleleAssessed().equals("C")) {
						eqtmz *= -1;
					}
					
					if (eqtmz < 0) {
						// flip all the zscores
						me.setZscore(me.getZscore() * -1);
						Double[] z = me.getDatasetZScores();
						for (int i = 0; i < z.length; i++) {
							if (z[i] != null) {
								z[i] *= -1;
							}
						}
					}
					
					String[] ds = me.getDatasets();
					for (int i = 0; i < ds.length; i++) {
						if (ds[i] != null && !ds[i].equals("-")) {
							ds[i] = "Meth-" + ds[i];
						}
					}
					
					outf.writeln(me.toString());
				}
				
				printed.add(eqtl);
				
			}
		}
		outf.close();
		
		
	}
	
	private void linkmePRSToEQTM(String meprs, String eqtm, String out) throws IOException {
		QTLTextFile f = new QTLTextFile(eqtm, QTLTextFile.R);
		ArrayList<EQTL> meqtl = f.readList();
		f.close();
		
		HashMap<String, String> cgToGen = new HashMap<String, String>();
		for (EQTL e : meqtl) {
			cgToGen.put(e.getRsName(), e.getProbe());
		}
		
		QTLTextFile f2 = new QTLTextFile(meprs, QTLTextFile.R);
		ArrayList<EQTL> eprs = f2.readList();
		f.close();
		
		TextFile outf = new TextFile(out, TextFile.W);
		outf.writeln("Trait\tCG\tCGChr\tCGChrPos\teQTMGene");
		for (EQTL e : eprs) {
			String cg = e.getProbe();
			String gene = cgToGen.get(cg);
			if (gene != null) {
				outf.writeln(e.getRsName() + "\t" + e.getProbe() + "\t" + e.getProbeChr() + "\t" + e.getProbeChrPos() + "\t" + gene);
			}
		}
		outf.close();
	}
	
	private void compareTP(String eqtmfile, String meqtlfile, String tpAndTn, String out) throws IOException {
		ArrayList<EQTL> tp = new ArrayList<EQTL>();
		ArrayList<EQTL> tn = new ArrayList<EQTL>();
		
		
		TextFile tf = new TextFile(tpAndTn, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		
		while (elems != null) {
			String tps = elems[0];
			String snp = elems[1];
			String gene = elems[2];
			
			EQTL e = new EQTL();
			e.setRsName(snp);
			e.setProbe(gene);
			if (tps.equals("TP")) {
				tp.add(e);
			} else {
				tn.add(e);
			}
			
			
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		System.out.println(tp.size() + " TP");
		System.out.println(tn.size() + " TN");
		
		QTLTextFile f = new QTLTextFile(eqtmfile, QTLTextFile.R);
		ArrayList<EQTL> eqtm = f.readList();
		f.close();
		System.out.println(eqtm.size() + " eQTM ");
		
		f = new QTLTextFile(meqtlfile, QTLTextFile.R);
		ArrayList<EQTL> meqtl = f.readList();
		f.close();
		System.out.println(meqtl.size() + " meQTL");
		TextFile tpout = new TextFile(out + "combos-tp.txt", TextFile.W);
		TextFile tnout = new TextFile(out + "combos-tn.txt", TextFile.W);
		TextFile tpandtnout = new TextFile(out + "combos-tnplustp.txt", TextFile.W);
		
		for (int q = 0; q < 2; q++) {
			ArrayList<EQTL> eqtl;
			int ctrmeqtl = 0;
			int ctreqtm = 0;
			int ctrmeqtlsignificant = 0;
			if (q == 0) {
				eqtl = tp;
			} else {
				eqtl = tn;
			}
			
			int qtlctr = 0;
			
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
					ctreqtm++;
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
					if (q == 0) {
						tpout.writeln(e.getRsName() + "\t" + maxeqtm.getRsName());
					} else {
						tnout.writeln(e.getRsName() + "\t" + maxeqtm.getRsName());
					}
					tpandtnout.writeln(e.getRsName() + "\t" + maxeqtm.getRsName());
				}
				
				
				// prepare output
				if (maxmeqtl != null) {
					ctrmeqtl++;
					if (maxmeqtl.getFDR() < 0.05) {
						ctrmeqtlsignificant++;
					}
					
				}
				
				qtlctr++;
				System.out.print("\r" + qtlctr + "/" + eqtl.size());
				
			}
			
			System.out.println();
			if (q == 0) {
				System.out.println("TP: " + tp.size() + "\teqtm: " + ctreqtm + "\tmeqtl: " + ctrmeqtl + "\tsignificant " + ctrmeqtlsignificant);
			} else {
				System.out.println("TN: " + tp.size() + "\teqtm: " + ctreqtm + "\tmeqtl: " + ctrmeqtl + "\tsignificant " + ctrmeqtlsignificant);
			}
			
		}
		
		tpandtnout.close();
		tnout.close();
		tpout.close();
		
		
	}
	
	private void compare2(String eqtlfile, String eqtmfile, String meqtlfile, String tpAndTn, String compout) throws Exception {
		
		ArrayList<EQTL> eqtl = null;
		QTLTextFile f = null;
		if (eqtlfile.endsWith("eQTLsWithIdenticalDirecton.txt.gz")) {
			TextFile tf = new TextFile(eqtlfile, TextFile.R);
			String[] elems = tf.readLineElems(TextFile.tab);
			eqtl = new ArrayList<>();
			while (elems != null) {
				String snp = elems[0];
				String gene = elems[1];
				String alleles = elems[2];
				String assessed = elems[3];
				String zscorestr = elems[4];
				Double z = Double.parseDouble(zscorestr);
				
				EQTL e = new EQTL();
				
				e.setProbe(gene);
				e.setRsName(snp);
				
				e.setAlleles(alleles);
				e.setAlleleAssessed(assessed);
				e.setZscore(z);
				eqtl.add(e);
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
		} else {
			f = new QTLTextFile(eqtlfile, QTLTextFile.R);
			eqtl = f.readList();
			f.close();
		}
		
		
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
		zp.init(2, new String[]{"eQTL", "meQTL"}, false, plotname);
		
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
				double eqtmz = maxeqtm.getZscore();
				if (!maxeqtm.getAlleleAssessed().equals("C")) {
					eqtmz *= -1;
				}
				if (eqtmz < 0) {
					zscoreflipped2 *= -1;
				}

//				zp.draw(e.getZscore(), eqtmz, 0, 1);
				zp.draw(e.getZscore(), zscoreflipped2, 0, 1);
//				zp.draw(eqtmz, zscoreflipped2, 1, 2);
				
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
