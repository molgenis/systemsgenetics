package nl.umcg.westrah.binarymetaanalyzer.westrah.binarymetaanalyzer.posthoc;

import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import javax.xml.soap.Text;
import java.io.IOException;
import java.util.HashSet;
import java.util.TreeSet;

public class QTLFileFilter {
	
	public static void main(String[] args) {
		
		String str = "a;b;c;d;-;f;g;h;-;j";
		int[] occ = Strings.countOccurrences(str, Strings.semicolon, "-");
		System.out.println(occ[0]);
		System.out.println(occ[1]);
		
		String[] meh = Strings.subsplit(str, Strings.semicolon, 9, 10);
		for (String s : meh) {
			System.out.println(s);
		}
		
	}
	
	public void filter(String snp, String gene, String snpgene, String in, int minimalNrOfCohorts, String out) throws IOException {
		
		System.out.println("QTL File Filter");
		System.out.println("Input: " + in);
		System.out.println("Output: " + out);
		System.out.println("SNPs: " + snp);
		System.out.println("Genes: " + gene);
		System.out.println("SNP+Genes: " + snpgene);
		HashSet<String> snps = null;
		HashSet<String> genes = null;
		HashSet<Pair<String, String>> snpgenes = null;
		
		if (snp != null) {
			snps = new HashSet<>();
			TextFile tf = new TextFile(snp, TextFile.R);
			snps.addAll(tf.readAsArrayList());
			tf.close();
			System.out.println(snps.size() + " snps in file");
		}
		
		if (genes != null) {
			genes = new HashSet<>();
			TextFile tf = new TextFile(gene, TextFile.R);
			genes.addAll(tf.readAsArrayList());
			tf.close();
			System.out.println(genes.size() + " genes in file");
		}
		
		if (snpgene != null) {
			
			snps = new HashSet<>();
			genes = new HashSet<>();
			snpgenes = new HashSet<>();
			
			TextFile tf = new TextFile(snpgene, TextFile.R);
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				if (elems.length >= 2) {
					snps.add(elems[0]);
					genes.add(elems[1]);
					snpgenes.add(new Pair<String, String>(elems[0], elems[1]));
				}
				elems = tf.readLineElems(TextFile.tab);
			}
			
			tf.close();
			
			System.out.println(snps.size() + " snps to filter for");
			System.out.println(genes.size() + " genes to filter for");
			System.out.println(snpgenes.size() + " SNP/Gene combinations");
		}
		
		TextFile tf = new TextFile(in, TextFile.R, 16 * 1048576);
		TextFile tfo = new TextFile(out, TextFile.W, 16 * 1048576);
		String header = tf.readLine();
		String[] headerElems = header.split("\t");
		boolean reducedFileFormat = false;
		if (headerElems.length <= 8) {
			reducedFileFormat = true;
			System.out.println("WARNING: Reduced file format detected.\nWill not be able to perform filtering for number of cohorts.\nFor #cohort filtering, use SNP/Gene filter instead");
		}
		tfo.writeln(header);
		
		String ln = tf.readLine();
		int ctr = 0;
		HashSet<String> genesFound = new HashSet<>();
		HashSet<String> snpsFound = new HashSet<>();
		
		
		int written = 0;
		while (ln != null) {
			String substr = null;
			if (ln.length() > 1000) {
				substr = ln.substring(0, 1000);
			} else {
				substr = ln;
			}
			
			int snpcol = 1;
			int genecol = 4;
			
			if (reducedFileFormat) {
				snpcol = 1;
				genecol = 2;
			}
			
			boolean write = true;
			String[] elems = substr.split("\t");
			if (snps != null && !snps.contains(elems[snpcol])) {
				write = false;
			} else {
				snpsFound.add(elems[snpcol]);
			}
			if (genes != null && !genes.contains(elems[genecol])) {
				write = false;
			} else {
				genesFound.add(elems[genecol]);
			}
			if (snpgenes != null && write) {
				Pair<String, String> s = new Pair<String, String>(elems[snpcol], elems[genecol]);
				if (!snpgenes.contains(s)) {
					write = false;
				}
			}
			
			if (!reducedFileFormat && minimalNrOfCohorts > 1) {
				// check number of cohorts
				String[] lnelems = Strings.subsplit(ln, Strings.tab, 11, 12);
				String datasets = lnelems[0];
				int[] occ = Strings.countOccurrences(datasets, Strings.semicolon, "-");
				int nrmissing = occ[1];
				int total = occ[0];
				if (total - nrmissing < minimalNrOfCohorts) {
					write = false;
				}
			}
			
			if (write) {
				tfo.writeln(ln);
				written++;
				if (snpgenes != null && snpgenes.size() == written) {
					System.out.println("I seem to have found all requested SNP/Gene combinations. Stopping.");
					break;
				}
			}
			ln = tf.readLine();
			ctr++;
			if (ctr % 100000 == 0) {
				System.out.print(ctr + " lines processed. " + genesFound.size() + " genes and " + snpsFound.size() + " snps found sofar. " + written + " lines written\r");
			}
		}
		
		
		tf.close();
		tfo.close();
		System.out.println();
		System.out.println("Done.");
	}
	
}
