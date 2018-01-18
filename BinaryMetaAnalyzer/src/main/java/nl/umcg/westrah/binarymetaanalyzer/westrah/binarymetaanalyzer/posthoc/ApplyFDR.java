package nl.umcg.westrah.binarymetaanalyzer.westrah.binarymetaanalyzer.posthoc;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;

public class ApplyFDR {
	
	
	public void apply(String in, String out, double threshold, boolean onlywritesignificant, String fdrFile) throws IOException {
		
		
		// read reference
		System.out.println("In:\t" + in);
		System.out.println("Out:\t" + out);
		System.out.println("FDR Reference:\t" + fdrFile);
		TextFile ref = new TextFile(fdrFile, TextFile.R);
		String[] refheader = ref.readLineElems(Strings.tab); // header
		int reffdr = -1;
		int refpval = -1;
		for (int i = 0; i < refheader.length; i++) {
			if (refheader[i].toLowerCase().equals("pvalue")) {
				refpval = i;
			}
			if (refheader[i].toLowerCase().equals("fdr")) {
				reffdr = i;
			}
		}
		
		System.out.println("Pval col:\t" + refpval);
		System.out.println("FDR col:\t" + reffdr);
		
		String ln = ref.readLine();
		LinkedList<Double> pvals = new LinkedList<Double>();
		LinkedList<Double> fdrs = new LinkedList<Double>();
		
		int lnctr = 1;
		while (ln != null) {
			String[] pval = Strings.subsplit(ln, Strings.tab, refpval, refpval + 1);
			String[] fdr = Strings.subsplit(ln, Strings.tab, reffdr, reffdr + 1);
			
			if (fdr == null || fdr.length == 0) {
				System.out.println("No FDR column found for line: " + lnctr);
				System.out.println(pval[0]);
				System.out.println(reffdr);
				System.out.println(ln);
				System.out.println(Strings.countSeparatorOccurrences(ln, Strings.tab));
				System.exit(-1);
			}
			
			
			
			double fdrd = Double.parseDouble(fdr[0]);
			if (onlywritesignificant && fdrd > threshold) {
				break;
			}
			pvals.add(Double.parseDouble(pval[0]));
			fdrs.add(fdrd);
			
			
			lnctr++;
			ln = ref.readLine();
		}
		ref.close();
		
		ref.close();
		TextFile tf = new TextFile(in, TextFile.R, 1048576);
		String[] header = tf.readLineElems(TextFile.tab);
		boolean hasFDRCol = false;
		int fdrcol = -1;
		for (int i = 0; i < header.length; i++) {
			if (header[i].equals("FDR")) {
				fdrcol = i;
				hasFDRCol = true;
			}
		}
		
		
		TextFile tfo = null;
		if (!onlywritesignificant) {
			tfo = new TextFile(out + ".txt.gz", TextFile.W, 1048576);
		}
		
		TextFile tfosig = new TextFile(out + "-Significant-" + threshold + ".txt.gz", TextFile.W);
		if (hasFDRCol) {
			if (!onlywritesignificant) {
				tfo.writeln(Strings.concat(header, Strings.tab));
			}
			tfosig.writeln(Strings.concat(header, Strings.tab));
		} else {
			if (!onlywritesignificant) {
				tfo.writeln(Strings.concat(header, Strings.tab) + "\tFDR");
			}
			tfosig.writeln(Strings.concat(header, Strings.tab) + "\tFDR");
		}
		
		
		double currentPVal = pvals.get(0);
		int currentPvalIndex = 0;
		double currentFDR = fdrs.get(0);
		
		String inln = tf.readLine();
		int ctr = 0;
		while (inln != null) {
			String[] pvale = Strings.subsplit(inln, Strings.tab, 0, 1);
			Double pval = Double.parseDouble(pvale[0]);
			
			// check the pvalue against the next fdr value threshold
			while (pval > currentPVal) {
				if (currentPvalIndex < pvals.size()) {
					currentPVal = pvals.get(currentPvalIndex);
					currentFDR = fdrs.get(currentPvalIndex);
				} else {
					currentPVal = 2;
					currentFDR = 1d;
				}
				// get the next
				currentPvalIndex++;
			}
			
			if (hasFDRCol) {
				String[] elems = Strings.tab.split(inln);
				
				// replace fdr value
				elems[fdrcol] = "" + currentFDR;
			} else {
				inln += "\t" + currentFDR;
			}
			
			if (currentFDR < threshold) {
				tfosig.writeln(inln);
			}
			if (!onlywritesignificant) {
				tfo.writeln(inln);
			} else {
				if (currentFDR >= threshold) {
					break;
				}
			}
			ctr++;
			if (ctr % 10000 == 0) {
				System.out.print(ctr + " lines parsed.\r");
			}
			inln = tf.readLine();
		}
		
		System.out.println();
		System.out.println();
		
		tfosig.close();
		if (!onlywritesignificant) {
			tfo.close();
		}
		tf.close();
	}
	
}
