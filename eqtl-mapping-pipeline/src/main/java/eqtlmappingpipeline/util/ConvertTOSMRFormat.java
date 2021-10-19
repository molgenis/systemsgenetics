package eqtlmappingpipeline.util;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashMap;

public class ConvertTOSMRFormat {
	
	public static void main(String[] args) {
		String file = "D:\\snpqc\\2018-01-31-cis-eQTLsFDR-ProbeLevel-CohortInfoFixed-8-sorted.txt.gz";
		String file2 = "D:\\snpqc\\2018-01-31-cis-eQTLsFDR-ProbeLevel-CohortInfoFixed-8-sorted-smr.txt.gz";
		String qc = "D:\\snpqc\\2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added-rewrite.txt.gz";
		ConvertTOSMRFormat c = new ConvertTOSMRFormat();
		try {
			c.run(file, qc, file2);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void run(String efile, String mafqcfile, String outfile) throws IOException {
		
		// load maf per snp
		HashMap<String, Double> afApersnp = new HashMap<String, Double>();
		HashMap<String, Double> afBpersnp = new HashMap<String, Double>();
		HashMap<String, String> minorallelepersnp = new HashMap<String, String>();
		HashMap<String, String> allelespersnp = new HashMap<String, String>();
		TextFile tf1 = new TextFile(mafqcfile, TextFile.R);
		System.out.println("Loading MAF info from: " + mafqcfile);
		tf1.readLine();
		int lnctr2 = 0;
		String[] elems = tf1.readLineElems(TextFile.tab);
		while (elems != null) {

            /*
            0 SNP
            1 Alleles
            2 AlleleA
            3 NumA
            4 NumB
            5 AF(A)
            6 AF(B)
            7 MinorAllele
            8 MAF
             */
			
			String snp = elems[0];
			double maf = Double.parseDouble(elems[8]);
			if (!Double.isNaN(maf)) {
				String minor = elems[7];
				String alleles = elems[1];
				allelespersnp.put(snp, alleles);
				minorallelepersnp.put(snp, minor);
				
				double afA = Double.parseDouble(elems[5]);
				double afB = Double.parseDouble(elems[6]);
				
				afApersnp.put(snp, afA);
				afBpersnp.put(snp, afB);
			}
			elems = tf1.readLineElems(TextFile.tab);
			lnctr2++;
			if (lnctr2 % 1000000 == 0) {
				System.out.println(lnctr2 + " SNPs parsed sofar.");
			}
		}
		tf1.close();
		System.out.println("MAF info found for " + afApersnp.size() + " variants, out of " + lnctr2 + " total.");
		
		// SNP     Chr     BP      A1      A2      Freq    Probe   Probe_Chr       Probe_bp        Gene    Orientation     b       se      p
		TextFile tf = new TextFile(efile, TextFile.R);
		TextFile tfo = new TextFile(outfile, TextFile.W);
		
		String lnout = "SNP" +
				"\tChr" +
				"\tBP" +
				"\tA1" +
				"\tA2" +
				"\tFreq" +
				"\tProbe" +
				"\tProbeChr" +
				"\tProbe_bp" +
				"\tGene" +
				"\tOrientation" +
				"\tb" +
				"\tse" +
				"\tp";
		tfo.writeln(lnout);
		
		String[] header = tf.readLineElems(TextFile.tab);
		elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {

            /*
            0 PValue
            1 SNPName
            2 SNPChr
            3 SNPChrPos
            4 ProbeName
            5 ProbeChr
            6 ProbeCenterChrPos
            7 CisTrans
            8 SNPType
            9 AlleleAssessed
            10 OverallZScore
            11 DatasetsWhereSNPProbePairIsAvailableAndPassesQC
            12 DatasetsZScores
            13 DatasetsNrSamples
            14 IncludedDatasetsMeanProbeExpression
            15 IncludedDatasetsProbeExpressionVariance
            16 HGNCName
            17 IncludedDatasetsCorrelationCoefficient
            18 Meta-Beta (SE)
            19 Beta (SE)
            20 FoldChange
            21 FDR

            0 "SNP" +
            1    "\tChr" +
            2    "\tBP" +
            3    "\tA1" +
            4    "\tA2" +
            5    "\tFreq" +
            6    "\tProbe" +
            7    "\tProbeChr" +
            8    "\tProbe_bp" +
            9    "\tGene" +
            10    "\tOrientation" +
            11    "\tb" +
            12    "\tse" +
            13    "\tp";
             */
			
			String snp = elems[1];
			
			String assessed = elems[9];
			
			Double refAfA = afApersnp.get(snp);
			
			if (refAfA == null) {
				System.err.println("Error: " + snp + " not in MAF file?");
			} else {
				
				String refalleles = allelespersnp.get(snp);
				String[] refallelelems = refalleles.split("/");
				String refAlleleA = refallelelems[0];
				String refAlleleB = refallelelems[1];
				
				Boolean flip = BaseAnnot.flipalleles(refalleles, refAlleleA, elems[8], assessed);
				
				if (flip != null) {
					String betastr = elems[18];
					String[] betastrelems = Strings.whitespace.split(betastr);
					Double beta = Double.parseDouble(betastrelems[0]);

					String se = betastrelems[1].replace(")", "");
					se = se.replace("(", "");
					
					if (flip) {
						beta = -beta;
					}
					
					String[] output = new String[14];
					output[0] = elems[1];
					output[1] = elems[2];
					output[2] = elems[3];
					output[3] = refAlleleA;
					output[4] = refAlleleB;
					output[5] = "" + refAfA;
					output[6] = elems[4];
					output[7] = elems[5];
					output[8] = elems[6];
					output[9] = elems[4];
					output[10] = "NA";
					output[11] = "" + beta;
					output[12] = se;
					output[13] = elems[0];
					tfo.writeln(Strings.concat(output, Strings.tab));
				}
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tfo.close();
		tf.close();
		
	}
}
