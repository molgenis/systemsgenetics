package nl.umcg.westrah.binarymetaanalyzer;

import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

public class MatrixConverterCis {
	private static final String NR_CALLED = "5075";
	private static final String MAF = "0.5";
	private static final String HWE = "1.0";
	private static final String CALL_RATE = "1.0";
	
	public void run(String realdatafile, String permutedfile, String out) throws IOException {
		// FROM: A1BG IL3872817 rs767575 -0.65 0.52 -1.06 0.32 0.81 -0.87 0.03 0.77 -0.33 0.16
		// TO: A1BG IL3872817 rs767575 C T -1.01
		
		
		// ColNames.txt.gz:
		// ENSG00000180016
		// ENSG00000213424_ENSG00000264058
		
		// RowNames
		
		class SNP {
			String name;
			String riskAllele;
			String nonRiskAllele;
			ArrayList<String> probes;
			ArrayList<Float> z = new ArrayList<>();
			
			@Override
			public boolean equals(Object o) {
				if (this == o) return true;
				if (o == null || getClass() != o.getClass()) return false;
				
				SNP snp = (SNP) o;
				
				return name != null ? name.equals(snp.name) : snp.name == null;
			}
			
			@Override
			public int hashCode() {
				return name != null ? name.hashCode() : 0;
			}
			
			public void clear() {
				probes = new ArrayList<>();
				z = new ArrayList<>();
			}
		}
		
		HashMap<String, SNP> qtls = new HashMap<String, SNP>();
		
		// start with building an index
		for (int perm = 0; perm < 11; perm++) {
			System.out.println("Processing permutation: " + perm);
			// clear SNP object Z-scores and probes
			{
				Set<String> snps = qtls.keySet();
				for (String key : snps) {
					qtls.get(key).clear();
				}
			}
			
			TextFile tf = null;
			if (perm == 0) {
				/* REAL data format
			Gene Transcript-clust risk non_risk Z-score
A1BG IL3872817 rs767575 C T -1.01
			 */
				tf = new TextFile(realdatafile, TextFile.R);
				String header = tf.readLine(); // we don't need this
				
			} else {
				/*
			Permutation format:
			NO HEADER, NO alleles
			A1BG IL3872817 rs767575 -0.65 0.52 -1.06 0.32 0.81 -0.87 0.03 0.77 -0.33 0.16
			 */
				tf = new TextFile(permutedfile, TextFile.R);
			}
			String[] elems = tf.readLineElems(Strings.whitespace);
			int ctr = 0;
			while (elems != null) {
				String gene = elems[0];
				String probe = elems[1];
				String combo = Strings.cache(gene + "_" + probe);
				String snpname = Strings.cache(elems[2]);
				if (perm == 0) {
					String riskallele = Strings.cache(elems[3]);
					String nonriskallele = Strings.cache(elems[4]);
					Float z = Float.parseFloat(elems[5]);
					SNP snpobj = qtls.get(snpname);
					if (snpobj == null) {
						snpobj = new SNP();
						snpobj.riskAllele = riskallele;
						snpobj.nonRiskAllele = nonriskallele;
						snpobj.name = snpname;
						snpobj.z = new ArrayList<>();
						snpobj.probes = new ArrayList<>();
					}
					
					snpobj.z.add(z);
					snpobj.probes.add(combo);
					qtls.put(snpname, snpobj);
					
				} else {
					int permzcol = 3 + (perm - 1);
					SNP snpobj = qtls.get(snpname);
					if (snpobj == null) {
						System.out.println("Permutation has SNP that is not present in non permuted data");
						System.out.println("Rs Id: " + snpname);
						System.exit(-1);
					}
					Float z = Float.parseFloat(elems[permzcol]);
					snpobj.z.add(z);
					snpobj.probes.add(combo);
				}
				
				elems = tf.readLineElems(Strings.whitespace);
				ctr++;
				if (ctr % 10000 == 0) {
					System.out.print(ctr + "\tlines parsed.\r");
				}
				
			}
			tf.close();
			System.out.println();
			System.out.println("Done reading");
			System.out.println();
			System.out.println(qtls.size() + " unique SNPs");
			int nrQTL = 0;
			for (String s : qtls.keySet()) {
				nrQTL += qtls.get(s).z.size();
			}
			System.out.println(nrQTL + " QTL");
			
			
			// we have all the data we need
			// write probe list first (because easy ;)
			HashSet<String> probeList = new HashSet<String>();
			ArrayList<SNP> qtlArr = new ArrayList<>();
			Set<String> snps = qtls.keySet();
			for (String key : snps) {
				probeList.addAll(qtls.get(key).probes);
				qtlArr.add(qtls.get(key));
			}
			
			String filenameCols = "Dataset-ColNames.txt.gz";
			String filenameRows = "Dataset-RowNames.txt.gz";
			String filenameDat = "Dataset.dat";
			if (perm > 0) {
				filenameCols = "Dataset-PermutationRound-" + perm + "-ColNames.txt.gz";
				filenameRows = "Dataset-PermutationRound-" + perm + "-RowNames.txt.gz";
				filenameDat = "Dataset-PermutationRound-" + perm + ".dat";
			}
			
			TextFile out1 = new TextFile(out + filenameCols, TextFile.W);
			System.out.println("Writing: " + out + filenameCols);
			for (String probe : probeList) {
				out1.writeln(probe);
			}
			out1.close();
			System.out.println("Done");
			
			// now comes the trickier part
			BinaryFile bf = new BinaryFile(out + filenameDat, BinaryFile.W, 1048576);
			bf.writeInt(1); // this will be a cis- dataset
			TextFile out2 = new TextFile(out + filenameRows, TextFile.W, 1048576);
			
			String header = "SNP\tAlleles\tMinorAllele\tAlleleAssessed\tNrCalled\tMaf\tHWE\tCallRate";
			out2.writeln(header);
			System.out.println("Writing: " + out + filenameRows);
			System.out.println("Writing: " + out + filenameDat);
			int wctr = 0;
			for (SNP s : qtlArr) {
				// write the Z-scores (easy part)
				for (Float z : s.z) {
					bf.writeFloat(z);
				}
				
				// write the SNP index
				// SNP     Alleles MinorAllele     AlleleAssessed  NrCalled        Maf     HWE     CallRate
				// rs78311930      C/T     T       T       2767    0.07878568847126854     0.8956349159520639      1.0     4       ENSG00000151093 ENSG00000151092 ENSG00000077097 ENSG00000077092
				String snpstr = s.name
						+ "\t" + s.nonRiskAllele + "/" + s.riskAllele
						+ "\t" + s.riskAllele
						+ "\t" + NR_CALLED
						+ "\t" + MAF
						+ "\t" + HWE
						+ "\t" + CALL_RATE
						+ "\t" + s.probes.size()
						+ "\t" + Strings.concat(s.probes, Strings.tab);
				
				out2.writeln(snpstr);
				wctr++;
				if (wctr % 10000 == 0) {
					System.out.print(wctr + " snps written\r");
				}
			}
			System.out.println();
			System.out.println("Done.");
			bf.close();
			out2.close();
		}
		
		
	}
}
