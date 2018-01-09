package nl.umcg.westrah.binarymetaanalyzer.westrah.binarymetaanalyzer.posthoc;

import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;

import javax.xml.soap.Text;
import java.io.IOException;
import java.util.HashSet;
import java.util.TreeSet;

public class QTLFileFilter {
	
	public void filter(String snp, String gene, String snpgene, String in, String out) throws IOException {
		
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
			
			TextFile tf = new TextFile(gene, TextFile.R);
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
		
		TextFile tf = new TextFile(in, TextFile.R, 1048576);
		TextFile tfo = new TextFile(out, TextFile.W, 1048576);
		tfo.writeln(tf.readLine());
		
		String ln = tf.readLine();
		int ctr = 0;
		while (ln != null) {
			String substr = null;
			if (ln.length() > 100) {
				substr = ln.substring(0, 100);
			} else {
				substr = ln;
			}
			boolean write = true;
			String[] elems = substr.split("\t");
			if (snps != null && !snps.contains(elems[1])) {
				write = false;
			}
			if (genes != null && !genes.contains(elems[4])) {
				write = false;
			}
			if (snpgenes != null && write) {
				Pair<String, String> s = new Pair<String, String>(elems[1], elems[4]);
				if (!snpgenes.contains(s)) {
					write = false;
				}
			}
			if (write) {
				tfo.writeln(ln);
			}
			ln = tf.readLine();
			ctr++;
			if (ctr % 10000 == 0) {
				System.out.print(ctr + " lines processed\r");
			}
		}
		
		
		tf.close();
		tfo.close();
		
		System.out.println("Done.");
	}
	
	public class SNPGenePairSet extends TreeSet<Pair<String, String>> {
	
	}
	
}
