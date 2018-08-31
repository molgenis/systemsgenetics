/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper.converters;

import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.io.trityper.util.ChrAnnotation;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;

/**
 * @author harmjan
 */
public class TriTyperToVCF {
	
	boolean hasDosage = false;
	
	public void convert(String in, String out, String snpSelection) throws IOException {
		
		System.out.println("Reading: " + in);
		if (!out.endsWith("/")) {
			out += "/";
		}
		Gpio.createDir(out);
		System.out.println("Writing: " + out + "export.vcf.gz");
		TriTyperGenotypeData d = new TriTyperGenotypeData();
		d.load(in);
		SNPLoader loader = d.createSNPLoader(10000);
		
		hasDosage = loader.hasDosageInformation();
		// make a list of included individuals..
		
		String[] individuals = d.getIndividuals();
		ArrayList<String> indsIncluded = new ArrayList<String>();
		HashMap<String, Integer> indToId = new HashMap<String, Integer>();
		
		for (int i = 0; i < individuals.length; i++) {
			if (d.getIsIncluded()[i] != null && d.getIsIncluded()[i]) {
				indsIncluded.add(individuals[i]);
				indToId.put(individuals[i], i);
				
			}
		}
		
		System.out.println("Exporting data for " + indsIncluded.size() + " individuals");
		HashSet<String> selectedSNPs = null;
		if (snpSelection != null) {
			selectedSNPs = new HashSet<String>();
			TextFile s = new TextFile(snpSelection, TextFile.R);
			selectedSNPs.addAll(s.readAsArrayList());
			s.close();
			System.out.println(selectedSNPs.size() + " snps selected from file: " + snpSelection);
		}
		
		
		String[] snps = d.getSNPs();
		int snpCtr = 0;
		for (String snp : snps) {
			if (snpSelection == null || selectedSNPs.contains(snp)) {
				snpCtr++;
			}
		}
		
		System.out.println("About to convert " + snpCtr + " SNPs to VCF");
		
		TextFile vcfOut = new TextFile(out + "export.vcf.gz", TextFile.W);
		writeheader(vcfOut, indsIncluded);
		ProgressBar pb = new ProgressBar(snpCtr, "Converting SNPs:");
		int nrexcluded = 0;
		int nrincluded = 0;
		int nrlowmaf = 0;
		int nrlowcr = 0;
		for (int s = 0; s < snps.length; s++) {
			if (snpSelection == null || selectedSNPs.contains(snps[s])) {
				SNP obj = d.getSNPObject(s);
				loader.loadGenotypes(obj);
				Byte chr = d.getChr(s);
				
				if (chr != null && chr > 0 && chr < 24) {
					if (obj.getCR() < 0.95) {
						nrlowcr++;
					}
					if (obj.getMAF() < 0.01) {
						nrlowmaf++;
					}
					if (obj.getCR() > 0.95 && obj.getMAF() > 0.01) {
						if (hasDosage) {
							loader.loadDosage(obj);
						}
						writeln(obj, d.getIsIncluded(), vcfOut);
						nrincluded++;
						obj.clearGenotypes();
					} else {
						nrexcluded++;
					}
				}
				pb.iterate();
			}
		}
		pb.close();
		System.out.println("Done. " + nrincluded + " included, " + nrexcluded + " excluded: " + nrlowcr + " low cr, " + nrlowmaf + " low maf");
		loader.close();
		vcfOut.close();
	}
	
	private void writeheader(TextFile out, ArrayList<String> indsIncluded) throws IOException {
		out.writeln("##fileformat=VCFv4.0");
		
		
		DateFormat dateFormat = new SimpleDateFormat("yyyyMMdd");
		Date date = new Date();
		out.writeln("##fileDate=" + dateFormat.format(date));
		
		out.writeln("##source=TriTyperToVCF");
		
		out.writeln("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">");
		out.writeln("##INFO=<ID=MA,Number=.,Type=String,Description=\"MinorAllele\">");
		out.writeln("##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">");
		out.writeln("##INFO=<ID=HW,Number=.,Type=Float,Description=\"Hardy-Weinberg P\">");
		out.writeln("##INFO=<ID=CR,Number=.,Type=Float,Description=\"Callrate\">");
		out.writeln("##FILTER=<ID=fqc,Description=\"HWE below 0.001 or maf below 0.05 or callrate below 0.95\">");
		out.writeln("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
		if (hasDosage) {
			out.writeln("##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Dosage\">");
		}
		
		// get all included individuals..
		String indList = Strings.concat(indsIncluded.toArray(new String[0]), Strings.tab);
		out.writeln("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + indList);
		
	}
	
	private void writeln(SNP snp, Boolean[] indIncluded, TextFile out) throws IOException {
		
		double maf = snp.getMAF();
		double hwe = snp.getHWEP();
		double cr = snp.getCR();
		String fqc = "PASS";
		if (maf < 0.05 || hwe < 0.001 || cr < 0.95) {
			fqc = "fqc";
		}
		String infoStr = "NS=" + snp.nrCalled
				+ ";MA=" + BaseAnnot.toString(snp.getMinorAllele())
				+ ";AF=" + maf
				+ ";HW=" + hwe
				+ ";CR=" + cr;
		
		String formatStr = "GT";
		if (hasDosage) {
			formatStr += ":DS";
		}
		byte allele1 = snp.getAlleles()[0];
		byte allele2 = snp.getAlleles()[1];
		
		String allele2Str = BaseAnnot.toString(allele2);
		//if (cr < 1 && allele2Str != null) {
		//    allele2Str += ",.";
		//} else
		
		if (allele2Str == null) {
			allele2Str = ".";
		}
		
		StringBuilder ln = new StringBuilder();
		ln.append(ChrAnnotation.parseByte(snp.getChr()))
				.append("\t").append(snp.getChrPos())
				.append("\t").append(snp.getName())
				.append("\t").append(BaseAnnot.toString(allele1))
				.append("\t").append(allele2Str)
				.append("\t").append(1)
				.append("\t").append(fqc)
				.append("\t").append(infoStr)
				.append("\t").append(formatStr);
		
		
		StringBuilder gtStr = new StringBuilder();
		byte[] alleles1 = snp.getAllele1();
		byte[] alleles2 = snp.getAllele2();
		
		
		double[] dosage = snp.getDosageValues();
		for (int i = 0; i < alleles1.length; i++) {
			if (indIncluded[i]) {
				
				byte indAllele1 = alleles1[i];
				byte indAllele2 = alleles2[i];
				
				char allelecode1 = '.';
				char allelecode2 = '.';
				if (indAllele1 == allele1) {
					allelecode1 = '0';
				} else if (indAllele1 == allele2) {
					allelecode1 = '1';
				}
				
				if (indAllele2 == allele1) {
					allelecode2 = '0';
				} else if (indAllele2 == allele2) {
					allelecode2 = '1';
				}
				
				short gt1 = alleles1[i];
				if (!hasDosage) {
					gtStr.append("\t").append(allelecode1).append("|").append(allelecode2);
				} else {
					
					gtStr.append("\t").append(allelecode1).append("|").append(allelecode2).append(":").append(dosage[i]);
				}
			}
		}
		
		
		ln.append(gtStr.toString());
		
		// now export all genotypes
		out.writeln(ln.toString());
	}
}
