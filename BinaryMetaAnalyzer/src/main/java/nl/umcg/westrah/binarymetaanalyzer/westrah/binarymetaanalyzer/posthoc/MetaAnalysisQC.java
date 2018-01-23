package nl.umcg.westrah.binarymetaanalyzer.westrah.binarymetaanalyzer.posthoc;

import org.apache.commons.math3.util.CombinatoricsUtils;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.stats.Heterogeneity;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;

public class MetaAnalysisQC {
	
	
	public static void main(String[] args) {
		// test this
		
		MetaAnalysisQC q = new MetaAnalysisQC();
		String efile = "C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\eqtl\\eQTLsFDR0.05-PrunedLevel_2CohortFilter_20170925.txt.gz";
		String outfile = efile + "-AllComboComparison.txt";
		
		
		String eqtlfiles = "C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\eqtl\\eQTLsFDR0.05-PrunedLevel_2CohortFilter_ParalogueFilter_20171204.txt.gz," +
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\eqtl\\stegle\\trans_eQTL_replication_CD4_Stegle_20171227-fdrfilter.txt.gz," +
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\eqtl\\stegle\\trans_eQTL_replication_CD8_Stegle_20171227-fdrfilter.txt.gz," +
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\eqtl\\stegle\\trans_eQTL_replication_CD14_Stegle_20171227-fdrfilter.txt.gz," +
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\eqtl\\stegle\\trans_eQTL_replication_CD15_Stegle_20171227-fdrfilter.txt.gz," +
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\eqtl\\stegle\\trans_eQTL_replication_CD19_Stegle_20171227-fdrfilter.txt.gz," +
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\eqtl\\stegle\\trans_eQTL_replication_Mac_Stegle_20171227-fdrfilter.txt.gz," +
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\eqtl\\stegle\\trans_eQTL_replication_PLT_Stegle_20171227-fdrfilter.txt.gz," +
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\eqtl\\2017-12-20-LCLMeta\\eQTLs-FDR0.05-fdrfilter.txt," +
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\eqtl\\methylationQTL-FDR0.05.txt.gz";
		
		
		String cohorttable = "C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\2017-12-27-CohortGroupsPlatforms.txt";
		String snpannotation = "";
		outfile = "C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\Heterogeneity\\2017-12-27-ZScoreComparisonTable.txt";
		boolean onlyconsidereffectsFromFirstFile = true;
		double threshold = 0.05;
		try {
			q.compareEffectSizes(eqtlfiles, cohorttable, outfile, null, threshold, false, onlyconsidereffectsFromFirstFile);
		} catch (IOException e) {
			e.printStackTrace();
		}

//		try {
////			q.effectSizeStability(efile, outfile);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
		
	}

//	public Pair<Double,Double> ZToBetaSE(double variance, double p, double q, int n, double z) {
//		double se = Math.sqrt(variance/(n*2pq));
//		double beta = se*z;
//		return new Pair<Double, Double>(beta,se);
//	}
	
	public void leaveOneOut(String eFile, String outfile) throws IOException {
		
		TextFile out = new TextFile(outfile, TextFile.W);
		
		String header = "SNP\tGene\tN\tiSquare\tP(iSquare)\tnumDatasets\tR\tRlo\tRavg\tRmedian\tRvar\tRhi";
		out.writeln(header);
		TextFile tf = new TextFile(eFile, TextFile.R);
		int nrLines = tf.countLines();
		System.out.println(nrLines + " in file: " + eFile);
		System.out.println("Output: " + outfile);
		tf.close();
		tf.open();
		tf.readLine();
		String[] data = tf.readLineElems(Strings.tab);
		ProgressBar pb = new ProgressBar(nrLines, "Performing leave one out meta-analysis...");
		while (data != null) {
			if (data.length > 13) {
				String snp = data[1];
				String gene = data[4];
				String alleles = data[8];
				String alleleAssessed = data[9];
				String[] cohorts = data[11].split(";");
				String[] cohortsZ = data[12].split(";");
				String[] cohortsNr = data[13].split(";");
				
				ArrayList<Double> z = new ArrayList<Double>(cohorts.length);
				ArrayList<Integer> n = new ArrayList<Integer>(cohorts.length);
				int nTotal = 0;
				for (int i = 0; i < cohortsZ.length; i++) {
					if (!cohorts[i].equals("-")) {
						z.add(Double.parseDouble(cohortsZ[i]));
						int nc = Integer.parseInt(cohortsNr[i]);
						n.add(nc);
						nTotal += nc;
					}
				}
				
				ArrayList<Double> rmetas = new ArrayList<Double>();
				for (int i = 0; i < z.size(); i++) {
					// select all datasets, except i;
					double[] q = new double[z.size() - 1];
					int[] r = new int[z.size() - 1];
					int ctr = 0;
					int metatotaln = 0;
					for (int j = 0; j < z.size(); j++) {
						if (j != i) {
							q[ctr] = z.get(j);
							r[ctr] = n.get(j);
							metatotaln += r[ctr];
							ctr++;
						}
					}
					
					// meta
					double zmeta = ZScores.getWeightedZ(q, r);
					double rmeta = ZScores.zToR(zmeta, metatotaln);
					rmetas.add(rmeta);
					
				}
				
				double[] zarr = Primitives.toPrimitiveArr(z);
				int[] narr = Primitives.toPrimitiveArr(n);
				int nmeta = 0;
				for (int q : narr) {
					nmeta += q;
				}
				
				double zmeta = ZScores.getWeightedZ(zarr, narr);
				double R = ZScores.zToR(zmeta, nmeta);
				// avg, etc
				// "SNP\tGene\tN\tnumDatasets\tR\tRlo\tRavg\tRmedian\tRvar\tRhi";
				
				double[] rmetaArr = Primitives.toPrimitiveArr(rmetas);
				double rlo = JSci.maths.ArrayMath.min(rmetaArr);
				double rhi = JSci.maths.ArrayMath.max(rmetaArr);
				double ravg = JSci.maths.ArrayMath.mean(rmetaArr);
				double rmed = JSci.maths.ArrayMath.median(rmetaArr);
				double rvar = JSci.maths.ArrayMath.variance(rmetaArr);
				int nDatasets = z.size();
				Triple<Double, Double, Integer> isqvals = Heterogeneity.getISq(z, n);
				double isquare = isqvals.getLeft();
				double isquarep = isqvals.getMiddle();
				String outln = snp + "\t" + gene
						+ "\t" + nmeta
						+ "\t" + isquare
						+ "\t" + isquarep
						+ "\t" + nDatasets
						+ "\t" + R
						+ "\t" + rlo
						+ "\t" + ravg
						+ "\t" + rmed
						+ "\t" + rvar
						+ "\t" + rhi;
				
				out.writeln(outln);
			}
			pb.iterate();
			data = tf.readLineElems(Strings.tab);
		}
		pb.close();
		out.close();
		tf.close();
		
	}
	
	public void effectSizeStability(String eQTLFileStr, String outFileStr) throws IOException {
		
		
		TextFile tf = new TextFile(eQTLFileStr, TextFile.R);
		int nrLines = tf.countLines();
		System.out.println(nrLines + " lines in eQTL file: " + eQTLFileStr);
		tf.close();
		tf.open();
		String header = tf.readLine();
		String[] data = tf.readLineElems(TextFile.tab); // header
		String headerOut = "snp\tgene\tmetaz\tmetar\tn\tminR\tavgR\tmaxR\tvarR\tNrCombos";
		TextFile out = new TextFile(outFileStr, TextFile.W);
		out.writeln(headerOut);
		int qtlctr = 0;
		ProgressBar pb = new ProgressBar(nrLines);
		while (data != null) {
			if (data.length > 13) {
				String snp = data[1];
				String gene = data[4];
				String alleles = data[8];
				String alleleAssessed = data[9];
				String[] cohorts = data[11].split(";");
				String[] cohortsZ = data[12].split(";");
				String[] cohortsNr = data[13].split(";");
				
				ArrayList<Double> z = new ArrayList<Double>(cohorts.length);
				ArrayList<Integer> n = new ArrayList<Integer>(cohorts.length);
				int nTotal = 0;
				for (int i = 0; i < cohortsZ.length; i++) {
					if (!cohorts[i].equals("-")) {
						z.add(Double.parseDouble(cohortsZ[i]));
						int nc = Integer.parseInt(cohortsNr[i]);
						n.add(nc);
						nTotal += nc;
					}
				}
				
				int nrDatasetsTotal = z.size();
				ArrayList<Double> r = new ArrayList<Double>();
				int comboctr = 0;
				for (int nrDatasetsInCombo = 2; nrDatasetsInCombo < nrDatasetsTotal + 1; nrDatasetsInCombo++) {
					Iterator<int[]> combos = CombinatoricsUtils.combinationsIterator(nrDatasetsTotal, nrDatasetsInCombo);
					
					while (combos.hasNext()) {
						int[] combo = combos.next();
						ArrayList<Double> selectZ = new ArrayList<Double>(combo.length);
						ArrayList<Integer> selectN = new ArrayList<Integer>(combo.length);
						
						int sumN = 0;
						for (int i : combo) {
							selectZ.add(z.get(i));
							selectN.add(n.get(i));
							sumN += n.get(i);
						}
						
						if (combo.length == selectN.size()) {
							// meta-analyze
							int[] n2 = Primitives.toPrimitiveArr(selectN);
							double[] z2 = Primitives.toPrimitiveArr(selectZ);
							double zm = ZScores.getWeightedZ(z2, n2);
							
							double rm = ZScores.zToR(zm, sumN);
							r.add(rm);
						}
						comboctr++;
						if (comboctr % 10000000 == 0) {
							System.out.print("\r" + qtlctr + "\t" + nrDatasetsInCombo + "/" + nrDatasetsTotal + "\t" + comboctr + "\t" + r.size());
						}
					}
				}
				System.out.println();
				
				double[] ra = Primitives.toPrimitiveArr(r);
				double[] za = Primitives.toPrimitiveArr(z);
				int[] na = Primitives.toPrimitiveArr(n);
				double min = JSci.maths.ArrayMath.min(ra);
				double max = JSci.maths.ArrayMath.max(ra);
				double avg = JSci.maths.ArrayMath.mean(ra);
				double var = JSci.maths.ArrayMath.variance(ra);
				double metaZ = ZScores.getWeightedZ(za, na);
				double metaR = ZScores.zToR(metaZ, nTotal);
				
				
				String outln = snp + "\t" + gene + "\t" + metaZ + "\t" + metaR + "\t" + nTotal + "\t" + min + "\t" + avg + "\t" + max + "\t" + var + "\t" + ra.length;
				out.writeln(outln);
				out.flush();
				
				qtlctr++;
				pb.set(qtlctr);
			}
			data = tf.readLineElems(TextFile.tab); // header
		}
		tf.close();
		out.close();
		pb.close();
		
	}
	
	/**
	 * @param eQTLFileStr         a string pointing to a (comma separated list of) eQTLFile(s)
	 * @param groupDefinitionFile a string pointing to a file containing two columns: cohortname\tgroupName
	 * @param outfileLoc          a string pointing to a location where I can save my output
	 * @param snpannotationfile   a string pointing to a file containing some annotation for the snps
	 * @throws IOException
	 */
	public void compareEffectSizes(String eQTLFileStr,
								   String groupDefinitionFile,
								   String outfileLoc,
								   String snpannotationfile,
								   double threshold,
								   boolean nonan,
								   boolean onlyConsiderEffectsFromFirstFile) throws IOException {
		
		// read the group definitions
		
		HashMap<String, ArrayList<Integer>> cohortToGroupIndex = new HashMap<String, ArrayList<Integer>>(); // point a cohort to a certain groupId
		HashMap<String, Integer> groupIndex = new HashMap<String, Integer>(); // point a group name to a certain groupId
		ArrayList<String> groupNames = new ArrayList<String>(); // a list of groupnames
		
		TextFile tf = new TextFile(groupDefinitionFile, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		int ctr = 0;
		while (elems != null) {
			if (elems.length > 1) {
				String cohort = elems[0];
				String group = elems[1];
				
				if (!groupIndex.containsKey(group)) {
					groupIndex.put(group, ctr);
					groupNames.add(group);
					System.out.println("New Group: " + group + ". ID: " + ctr);
					ctr++;
				}
				Integer id = groupIndex.get(group);
				
				ArrayList<Integer> groupIds = cohortToGroupIndex.get(cohort);
				if (groupIds == null) {
					groupIds = new ArrayList<>();
				}
				groupIds.add(id);
				
				cohortToGroupIndex.put(cohort, groupIds);
				System.out.println("Adding cohort " + cohort + " to group " + group + " with id " + groupIndex.get(group));
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		
		int nrgroups = groupIndex.size();
		int[] cohortspergroup = new int[nrgroups];
		HashMap<String, Integer> cohortIndexWithInGroup = new HashMap<String, Integer>();
		{
			int[] groupDsCtr = new int[nrgroups];
			for (String cohort : cohortToGroupIndex.keySet()) {
				ArrayList<Integer> gid = cohortToGroupIndex.get(cohort);
				for (int g : gid) {
					cohortIndexWithInGroup.put(cohort + "-" + g, groupDsCtr[g]);
					groupDsCtr[g]++;
					cohortspergroup[g]++;
				}
				
			}
		}
		System.out.println(nrgroups + " groups defined..");
		System.out.println(cohortToGroupIndex.size() + " cohorts specified.");
		
		// read in the QTL
		
		// isq per file
		// isq per group
		boolean debug = false;
		String[] eQTLFiles = eQTLFileStr.split(",");
		HashMap<QTLQC, QTLQC> qtlmap = new HashMap<QTLQC, QTLQC>(); // this feels like inception
		for (int q = 0; q < eQTLFiles.length; q++) {
			String eQTLFile = eQTLFiles[q];
			System.out.println("Parsing eQTL file: " + eQTLFile);
			TextFile in = new TextFile(eQTLFile, TextFile.R);
			String str = in.readLine(); // header
			int lnctr = 0;
			while ((str = in.readLine()) != null) {
				
				String[] data = str.split("\t");
				if (data.length > 13) {
					String snp = data[1];
					String snpChr = data[2];
					String snpChrPos = data[3];
					String gene = data[4];
					String geneChr = data[5];
					String geneChrPos = data[6];
					String alleles = data[8];
					String alleleAssessed = data[9];
					String[] cohorts = data[11].split(";");
					String[] cohortsZ = data[12].split(";");
					String[] cohortsNr = data[13].split(";");
					
					
					// see if we've already seen this snp..
					QTLQC qtl = new QTLQC(snp, gene, alleles, alleleAssessed, nrgroups, eQTLFiles.length, cohortspergroup);
					Boolean flipAlleles = null;
					
					boolean add = true;
					if (onlyConsiderEffectsFromFirstFile && q > 0) {
						if (!qtlmap.containsKey(qtl)) {
							add = false;
						}
					}
					
					if (qtlmap.containsKey(qtl)) {
						qtl = qtlmap.get(qtl);
						String allelesOrig = qtl.getAlleles();
						String allelesAssessedOrig = qtl.getAlleleAssessed();
						flipAlleles = BaseAnnot.flipalleles(allelesOrig, allelesAssessedOrig, alleles, alleleAssessed);
					} else {
						if (add) {
							qtlmap.put(qtl, qtl);
						}
						flipAlleles = false;
					}
					
					if (flipAlleles == null) {
						System.err.println("Warning:  " + snp + "/" + gene + " eQTL found alleles " + alleles + ", assessed " + alleleAssessed + " while expecting: " + qtl.getAlleles() + ", " + qtl.getAlleleAssessed());
					}
					
					
					if (debug && snp.equals("rs28498283") && gene.equals("ENSG00000085265")) {
						System.out.println(qtl);
//						debug = true;
						System.out.println("flipperdeflip " + flipAlleles);
					}
					
					qtl.setChr(snpChr, snpChrPos, geneChr, geneChrPos);
					ArrayList<Double> isqz = new ArrayList<Double>();
					ArrayList<Integer> isqn = new ArrayList<Integer>();
					int totaln = 0;
					for (int i = 0; i < cohorts.length; i++) {
						if (!cohorts[i].equals("-")) {
							Double z = Double.parseDouble(cohortsZ[i]);
							Integer n = Integer.parseInt(cohortsNr[i]);
							isqz.add(z);
							isqn.add(n);
							totaln += n;
						}
					}
					
					Triple<Double, Double, Integer> isq = Heterogeneity.getISq(isqz, isqn); // i,p
					double efilezmeta = ZScores.getWeightedZ(Primitives.toPrimitiveArr(isqz.toArray(new Double[0])),
							Primitives.toPrimitiveArr(isqn.toArray(new Integer[0])));
					double efilecor = ZScores.zToR(efilezmeta, isq.getRight());
					qtl.setQTLFileISq(q, isq.getMiddle(), isq.getLeft(), isq.getRight());
					qtl.setQTLFileCor(efilecor, q);

//					qtl.setAlleleflip(flipAlleles, eqtlfileId);
					
					
					for (int c = 0; c < cohorts.length; c++) {
						String cohortName = cohorts[c];
						ArrayList<Integer> groupIds = cohortToGroupIndex.get(cohortName);
						if (groupIds != null) {
							for (Integer groupId : groupIds) {
								if (debug) {
									System.out.println(q + "\t" + cohortName + "\t" + groupId);
								}
								if (groupId != null) {
									double z = Double.parseDouble(cohortsZ[c]);
									
									if (flipAlleles != null && flipAlleles) {
										z = -z;
									}
									
									int n = Integer.parseInt(cohortsNr[c]);
									qtl.addz(groupId, z, n);
									qtl.addn(groupId, n);
									Integer cohortIndexInGroup = cohortIndexWithInGroup.get(cohortName + "-" + groupId);
									qtl.setCohortZScore(groupId, cohortIndexInGroup, z, n);
									if (debug) {
										System.out.println(q + "\t" + cohortName + "\t" + groupId + "\t" + z + "\t" + n);
									}
								}
							}
						}
					}
					
				}
				lnctr++;
				if (lnctr % 1000 == 0) {
					System.out.print(lnctr + " lines parsed\r");
				}
			}
			in.close();
			System.out.println("Done. " + qtlmap.size() + " QTLs found sofar");
		}
		
		// read snp annotation, if any
		HashMap<String, String> hashSNPAnnotation = null;
		if (snpannotationfile != null) {
			hashSNPAnnotation = new HashMap<String, String>();
			System.out.println("Loading SNP annotation from: " + snpannotationfile);
			TextFile in = new TextFile(snpannotationfile, TextFile.R);
			String str = "";
			while ((str = in.readLine()) != null) {
				String[] data = str.split("\t");
				hashSNPAnnotation.put(data[0], str);
			}
			in.close();
		}
		
		
		String header = "SNP\tGene";
		for (int group1 = 0; group1 < groupIndex.size(); group1++) {
			String gname1 = groupNames.get(group1);
			header += "\tZScore-" + gname1 + "\tN-" + gname1;
		}
		
		TextFile tfout = new TextFile(outfileLoc, TextFile.W);
		tfout.writeln(header);
		// get all QTLs
		ArrayList<QTLQC> qtlList = new ArrayList<QTLQC>();
		qtlList.addAll(qtlmap.keySet());
		
		for (QTLQC e : qtlList) {
			String ln = e.getSnp() + "\t" + e.getGene();
			for (int group = 0; group < groupIndex.size(); group++) {
				Double z = e.getZ()[group];
				if (z == null) {
					ln += "\tNS\t" + e.getN()[group];
				} else {
					z /= Math.sqrt(e.getN()[group]);
					ln += "\t" + z + "\t" + e.getN()[group];
				}
				
			}
			tfout.writeln(ln);
		}
		
		tfout.close();

//		// construct a header
//		String header = "SNP\tSNPChr\tSNPChrPos\tGene\tGeneChr\tGeneChrPos";
//		for (int i = 0; i < eQTLFiles.length; i++) {
//			header += "\tCor-eQTLFile-" + i + "\tIsq-eQTLFile-" + i + "\tIsqP-eQTLFile-" + i + "\tIsqN-eQTLFile-" + i;
//		}
//		for (int group1 = 0; group1 < groupIndex.size(); group1++) {
//			String gname1 = groupNames.get(group1);
//			header += "\tCorr-" + gname1 + "\tN-" + gname1 + "\tIsq-" + gname1 + "\tIsqP-" + gname1 + "\tIsqN-" + gname1;
//		}
//		for (int group1 = 0; group1 < groupIndex.size(); group1++) {
//			String gname1 = groupNames.get(group1);
//			for (int group2 = group1 + 1; group2 < groupIndex.size(); group2++) {
//				String gname2 = groupNames.get(group2);
//				header += "\tZDiff-" + gname1 + "-" + gname2 + "\tPVal-" + gname1 + "-" + gname2;
//			}
//		}
//
//		header += "\tNrTests\tNrTestsSignificant";
//		if (hashSNPAnnotation != null) {
//			// format: rs7923609	10	65133822	A/G	G	Liver enzyme levels (alkaline phosphatase)	22001757
//
//			header += "\tRSId\tChr\tPos\tAlleles\tAssessed\tTraits\tPMIDs";
//		}
//
//		// prepare the outfile
//		TextFile tfout = new TextFile(outfileLoc, TextFile.W);
//		tfout.writeln(header);
//
//		// get all QTLs
//		ArrayList<QTLQC> qtlList = new ArrayList<QTLQC>();
//		qtlList.addAll(qtlmap.keySet());
//		if (debug) {
//			for (QTLQC e : qtlList) {
//				if (e.getSnp().equals("rs28498283") && e.getGene().equals("ENSG00000085265")) {
//					System.out.println(e);
//				}
//			}
//		}
//
//		// iterate
//		System.out.println("Comparing effect sizes for " + qtlList.size() + " QTL.");
//		int written = 0;
//		int[][] writtenArr = new int[nrgroups][nrgroups];
//		int[][] sigArr = new int[nrgroups][nrgroups];
//		int[] nrCompsNonSig = new int[qtlList.size()];
//
//		for (int q = 0; q < qtlList.size(); q++) {
//			double sumR21 = 0;
//			double sumR22 = 0;
//
//			QTLQC qtl = qtlList.get(q);
//			for (int groupId = 0; groupId < nrgroups; groupId++) {
//				qtl.dividezBySqrtN(groupId);
//			}
//
//			// check if qtl has any samples involved
//			int totalN = qtl.getTotalN();
//
//			// check if all the samples are from 1 cohort
////			int[] n = qtl.getN();
////			boolean allNFromSingleDS = qtl.getIsTotalNFromSingleGroup();
////
//
//			if (totalN > 0) {
////			for (int p = 0; p < nrTransPerSNP[snp]; p++) {
//				String lnOut = qtl.getSnp() + "\t" + qtl.getSnpChr() + "\t" + qtl.getSnpChrPos()
//						+ "\t" + qtl.getGene() + "\t" + qtl.getGeneChr() + "\t" + qtl.getSnpChrPos();
//
//				// get ISQ per eQTL file
//				Triple<double[], double[], int[]> efileIsq = qtl.getIsqEfile();
//				double[] isqEFile = efileIsq.getLeft();
//				double[] isqpEFile = efileIsq.getMiddle();
//				int[] isqnEFile = efileIsq.getRight();
//				double[] corEfile = qtl.getQTLFileCor();
//				for (int i = 0; i < isqEFile.length; i++) {
//					lnOut += "\t" + corEfile[i] + "\t" + isqEFile[i] + "\t" + isqpEFile[i] + "\t" + isqnEFile[i];
//				}
//
//				// add iSq per group
//				Triple<double[], double[], int[]> groupIsq = qtl.getGroupIsq();
//				double[] isqgroup = groupIsq.getLeft();
//				double[] isqpgroup = groupIsq.getMiddle();
//				int[] isqngroup = groupIsq.getRight();
//
//				double[] zs = new double[groupIndex.size()];
//				int[] ns = new int[groupIndex.size()];
//
//				for (int group1 = 0; group1 < nrgroups; group1++) {
//					double z1orig = qtl.getZ()[group1]; // snpZScores[group1][snp][p];
//					int n1 = qtl.getN()[group1]; // snpZScoresNr[group1][snp][p];
//					double corr1 = ZScores.zToR(z1orig, n1);
//					double z1 = zFromCorr(corr1);
//					zs[group1] = z1;
//					ns[group1] = n1;
//					lnOut += "\t" + corr1 + "\t" + n1 + "\t" + isqgroup[group1] + "\t" + isqpgroup[group1] + "\t" + isqngroup[group1];
//				}
//
//				// format :
//				// group1corr\tgroup1n\tgroup2corr\tgroup2n\tzdiff\tpval, etc
//				ArrayList<Boolean> isNaN = new ArrayList<Boolean>();
//				int[][] writtenArrTmp = new int[nrgroups][nrgroups];
//				int[][] sigArrTmp = new int[nrgroups][nrgroups];
//				int nrTests = 0;
//				int nrSignificantTests = 0;
//				for (int group1 = 0; group1 < nrgroups; group1++) {
//					int n1 = ns[group1];
//					double z1 = zs[group1];
//					for (int group2 = group1 + 1; group2 < nrgroups; group2++) {
//						int n2 = ns[group2];
//						double z2 = zs[group2];
//						double se = Math.sqrt((1 / (n1 - 3d)) + (1 / (n2 - 3d)));
//						double zDiff = (z1 - z2) / se;
//						double pVal = cern.jet.stat.Probability.normal(-Math.abs(zDiff)) * 2d;
//						boolean pisnan = Double.isNaN(pVal);
//						isNaN.add(pisnan);
//						if (!pisnan) {
//							nrTests++;
//							writtenArrTmp[group1][group2]++;
//							if (pVal < threshold) {
//								sigArrTmp[group1][group2]++;
//								nrSignificantTests++;
//							} else {
//								nrCompsNonSig[q]++;
//							}
//						}
//						lnOut += "\t" + zDiff + "\t" + pVal;
//					}
//				}
//
//				lnOut += "\t" + nrTests + "\t" + nrSignificantTests;
//				if (hashSNPAnnotation != null) {
//					lnOut += "\t" + hashSNPAnnotation.get(qtl.getSnp());
//				}
//
//				// check if all results are NaN;
//				int nanctr = 0;
//				for (Boolean b : isNaN) {
//					if (b) nanctr++;
//				}
////				if (nanctr != isNaN.size()) {
//
//				if (!nonan || nanctr == 0) {
//					// update ctrs
//					for (int group1 = 0; group1 < nrgroups; group1++) {
//						for (int group2 = group1 + 1; group2 < nrgroups; group2++) {
//							writtenArr[group1][group2] += writtenArrTmp[group1][group2];
//							sigArr[group1][group2] += sigArrTmp[group1][group2];
//						}
//					}
//
//					tfout.writeln(lnOut);
//					written++;
//				}
//
//			}
//
//		}
//		tfout.close();
//
//
//		System.out.println(written + " QTL written eventually compared.");
//
//		System.out.println();
//		System.out.println();
//		String header2 = "Group1\tGroup2\tNrQTL\tNrP<0.05\tPerc";
//		System.out.println(header2);
//
//		TextFile tf2 = new TextFile(outfileLoc + "_comparisons.txt", TextFile.W);
//		tf2.writeln(header2);
//		int nrcomps = 0;
//		for (int group1 = 0; group1 < groupIndex.size(); group1++) {
//			for (int group2 = group1 + 1; group2 < groupIndex.size(); group2++) {
//				String ln = groupNames.get(group1)
//						+ "\t" + groupNames.get(group2)
//						+ "\t" + writtenArr[group1][group2]
//						+ "\t" + sigArr[group1][group2]
//						+ "\t" + ((double) sigArr[group1][group2] / writtenArr[group1][group2]);
//				System.out.println(ln);
//				tf2.writeln(ln);
//				nrcomps++;
//			}
//		}
//		tf2.close();

//		int[] nrSig = new int[nrcomps];
//		for (int q = 0; q < nrCompsNonSig.length; q++) {
//			nrSig[nrCompsNonSig[q]]++;
//		}
//
//		System.out.println("NrTests\tSignificantly different in NrTests\t");
//		for (int q = 0; q < nrSig.length; q++) {
//			System.out.println(q + "\t" + nrSig[q]);
//		}
		
	}
	
	private double zFromCorr(double corr1) {
		double raplus = 1 * corr1 + 1;
		double raminus = 1 - corr1;
		double z = (Math.log(raplus) - Math.log(raminus)) / 2;
		return z;
	}
	
	
	private class QTLQC implements Comparable<QTLQC> {
		private final double[] isqFile;
		private final double[] isqPFile;
		String snp;
		String alleles;
		String alleleAssessed;
		String gene;
		
		boolean[] alleleflip;
		Double[] z;
		int[] n;
		
		double[][] groupzscores; // [group][datasets]
		int[][] groupsamplesizes;
		private int[] isqNFile;
		private double[] QTLFileCor;
		private String snpChr;
		private String snpChrPos;
		private String geneChr;
		private String geneChrPos;
		
		public QTLQC(String snp, String gene, String alleles, String alleleAssessed, int nrGroups, int nrEQTLFiles, int[] cohortsPerGroup) {
			this.z = new Double[nrGroups];
			this.n = new int[nrGroups];
			this.alleleflip = new boolean[nrGroups];
			this.snp = snp;
			this.gene = gene;
			this.alleles = alleles;
			this.alleleAssessed = alleleAssessed;
			this.isqFile = new double[nrEQTLFiles];
			this.isqPFile = new double[nrEQTLFiles];
			this.isqNFile = new int[nrEQTLFiles];
			this.QTLFileCor = new double[nrEQTLFiles];
			this.groupzscores = new double[nrGroups][];
			this.groupsamplesizes = new int[nrGroups][];
			
			for (int q = 0; q < nrGroups; q++) {
				int nrCohorts = cohortsPerGroup[q];
				groupzscores[q] = new double[nrCohorts];
				groupsamplesizes[q] = new int[nrCohorts];
				for (int z = 0; z < nrCohorts; z++) {
					groupzscores[q][z] = Double.NaN;
				}
			}
		}
		
		public void setQTLFileISq(int e, double p, double isq, int n) {
			this.isqFile[e] = isq;
			this.isqPFile[e] = p;
			this.isqNFile[e] = n;
		}
		
		public Triple<double[], double[], int[]> getIsqEfile() {
			return new Triple<double[], double[], int[]>(isqFile, isqPFile, isqNFile);
		}
		
		@Override
		public boolean equals(Object o) {
			if (this == o) return true;
			if (o == null || getClass() != o.getClass()) return false;
			
			QTLQC qtlqc = (QTLQC) o;
			
			if (!snp.equals(qtlqc.snp)) return false;
			return gene.equals(qtlqc.gene);
		}
		
		@Override
		public int hashCode() {
			int result = snp.hashCode();
			result = 31 * result + gene.hashCode();
			return result;
		}
		
		public String getSnp() {
			return snp;
		}
		
		public String getAlleles() {
			return alleles;
		}
		
		public String getAlleleAssessed() {
			return alleleAssessed;
		}
		
		public String getGene() {
			return gene;
		}
		
		public Double[] getZ() {
			return z;
		}
		
		public int[] getN() {
			return n;
		}
		
		public void setAlleleflip(boolean b, int group) {
			alleleflip[group] = b;
		}
		
		public void setSnp(String snp) {
			this.snp = snp;
		}
		
		public void setAlleles(String alleles) {
			this.alleles = alleles;
		}
		
		public void setAlleleAssessed(String alleleAssessed) {
			this.alleleAssessed = alleleAssessed;
		}
		
		public void setGene(String gene) {
			this.gene = gene;
		}
		
		public void setZ(Double[] z) {
			this.z = z;
		}
		
		public void setN(int[] n) {
			this.n = n;
		}
		
		public void addz(Integer groupId, double z) {
			this.z[groupId] += z;
		}
		
		public void addn(Integer groupId, int n) {
			this.n[groupId] += n;
		}
		
		public void dividezBySqrtN(int groupId) {
			this.z[groupId] /= Math.sqrt(n[groupId]);
		}
		
		@Override
		public int compareTo(QTLQC that) {
			if (this.equals(that)) {
				return 0;
			} else if (this.snp.equals(that.snp)) {
				// sort on probe
				return this.gene.compareTo(that.gene);
			}
			return 0;
		}
		
		public void addz(Integer groupId, double z, int n) {
			if (this.z[groupId] == null) {
				this.z[groupId] = 0d;
			}
			this.z[groupId] += (z * Math.sqrt(n));
		}
		
		public int getTotalN() {
			int sum = 0;
			for (int i : n) {
				sum += i;
			}
			return sum;
		}
		
		public boolean getIsTotalNFromSingleGroup() {
			int totalN = getTotalN();
			for (int i : n) {
				if (totalN == i) {
					return true;
				}
			}
			return false;
		}
		
		@Override
		public String toString() {
			return "QTLQC{" +
					"snp='" + snp + '\'' +
					", alleles='" + alleles + '\'' +
					", alleleAssessed='" + alleleAssessed + '\'' +
					", gene='" + gene + '\'' +
					", z=" + Arrays.toString(z) +
					", n=" + Arrays.toString(n) +
					'}';
		}
		
		public void setCohortZScore(Integer groupId, Integer cohortIndexInGroup, double z, int n) {
			this.groupzscores[groupId][cohortIndexInGroup] = z;
			this.groupsamplesizes[groupId][cohortIndexInGroup] = n;
		}
		
		public Triple<double[], double[], int[]> getGroupIsq() {
			
			double[] p = new double[groupsamplesizes.length];
			double[] i = new double[groupsamplesizes.length];
			int[] ntotal = new int[groupsamplesizes.length];
			for (int q = 0; q < groupsamplesizes.length; q++) {
				int nrNull = 0;
				for (int cohort = 0; cohort < groupsamplesizes[q].length; cohort++) {
					if (groupsamplesizes[q][cohort] == 0) nrNull++;
				}
				double[] z = new double[groupsamplesizes[q].length - nrNull];
				int[] n = new int[groupsamplesizes[q].length - nrNull];
				int ctr = 0;
				for (int cohort = 0; cohort < groupsamplesizes[q].length; cohort++) {
					if (groupsamplesizes[q][cohort] > 0) {
						z[ctr] = groupzscores[q][cohort];
						n[ctr] = groupsamplesizes[q][cohort];
						ctr++;
					}
				}
				Triple<Double, Double, Integer> isq = Heterogeneity.getISq(z, n);
				p[q] = isq.getMiddle();
				i[q] = isq.getLeft();
				ntotal[q] = isq.getRight();
			}
			
			
			return new Triple<double[], double[], int[]>(i, p, ntotal);
			
			
		}
		
		public void setQTLFileCor(double QTLFileCor, int e) {
			this.QTLFileCor[e] = QTLFileCor;
		}
		
		public double[] getQTLFileCor() {
			return QTLFileCor;
		}
		
		public void setChr(String snpChr, String snpChrPos, String geneChr, String geneChrPos) {
			this.snpChr = snpChr;
			this.snpChrPos = snpChrPos;
			this.geneChr = geneChr;
			this.geneChrPos = geneChrPos;
		}
		
		public String getSnpChr() {
			return snpChr;
		}
		
		public String getSnpChrPos() {
			return snpChrPos;
		}
		
		public String getGeneChr() {
			return geneChr;
		}
		
		public String getGeneChrPos() {
			return geneChrPos;
		}
	}
	
}
