package nl.umcg.westrah.binarymetaanalyzer.westrah.binarymetaanalyzer.posthoc;

import gnu.trove.map.hash.TDoubleDoubleHashMap;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

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
		TDoubleDoubleHashMap thresholds = new TDoubleDoubleHashMap();
		
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
			
			
			thresholds.put(Double.parseDouble(pval[0]), fdrd);
			if (fdrd >= 0.999) {
				System.out.println("Found FDR of " + fdrd + ": maximum reached. Breaking!");
				break;
			}
			
			lnctr++;
			if (lnctr % 100000 == 0) {
				System.out.print(lnctr + " lines parsed. " + thresholds.size() + " unique pvalues spotted. Current FDR: " + pval[0] + " --> " + fdrd + ".\r");
			}
			ln = ref.readLine();
		}
		ref.close();
		System.out.println();
		System.out.println(lnctr + " lines parsed in total.");
		System.out.println(thresholds.size() + " unique p-values");
		
		ref.close();
		TextFile tf = new TextFile(in, TextFile.R, 10 * 1048576);
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
		
		
		class DoublePair implements Comparable<DoublePair> {
			double a;
			double b;
			
			public DoublePair(double a, double b) {
				this.a = a;
				this.b = b;
			}
			
			@Override
			public boolean equals(Object o) {
				if (this == o) return true;
				if (o == null || getClass() != o.getClass()) return false;
				
				DoublePair that = (DoublePair) o;
				
				if (Double.compare(that.a, a) != 0) return false;
				return Double.compare(that.b, b) == 0;
			}
			
			@Override
			public int hashCode() {
				int result;
				long temp;
				temp = Double.doubleToLongBits(a);
				result = (int) (temp ^ (temp >>> 32));
				temp = Double.doubleToLongBits(b);
				result = 31 * result + (int) (temp ^ (temp >>> 32));
				return result;
			}
			
			
			@Override
			public int compareTo(DoublePair o) {
				if (this.equals(o)) {
					return 0;
				} else if (this.a > o.a) {
					return 1;
				} else {
					return -1;
				}
			}
		}
		
		// put pvals back to array
		double[] uniquepvals = thresholds.keys();
		ArrayList<DoublePair> pvalpairs = new ArrayList<>();
		for (double p : uniquepvals) {
			DoublePair z = new DoublePair(p, thresholds.get(p));
			pvalpairs.add(z);
		}
		Collections.sort(pvalpairs);
		
		System.out.println("Sort order: ");
		for (int q = 0; q < pvalpairs.size(); q++) {
			if (q < 10) {
				System.out.println(q + "\t" + pvalpairs.get(q).a);
			}
			if (q > 0 && pvalpairs.get(q).a < pvalpairs.get(q - 1).a) {
				System.out.println("Not correctly sorted!");
				System.out.println("Current pval: " + pvalpairs.get(q).a + " prev " + pvalpairs.get(q - 1).a);
				System.exit(-1);
			}
		}
		
		double currentPVal = pvalpairs.get(0).a;
		int currentPvalIndex = 0;
		double currentFDR = pvalpairs.get(0).b;
		
		System.out.println("Will replace FDRs now");
		String inln = tf.readLine();
		int ctr = 0;
		int significant = 0;
		while (inln != null) {
			String[] pvale = Strings.subsplit(inln, Strings.tab, 0, 1);
			Double pval = Double.parseDouble(pvale[0]);
			
			// check the pvalue against the next fdr value threshold
			while (pval > currentPVal) {
				if (currentPvalIndex < pvalpairs.size()) {
					currentPVal = pvalpairs.get(currentPvalIndex).a;
					currentFDR = pvalpairs.get(currentPvalIndex).b;
				} else {
					currentPVal = 2;
					currentFDR = 1d;
				}
				// get the next
				currentPvalIndex++;
			}
			
			if (hasFDRCol) {
				String[] elems = Strings.tab.split(inln);
				// check whether the columns is actually there
				// e.g. fdrcol = 21, but len = 20
				if (elems.length <= fdrcol) {
					inln += "\t" + currentFDR;
				} else {
					// replace fdr value
					elems[fdrcol] = "" + currentFDR;
					inln = Strings.concat(elems, Strings.tab);
				}
			} else {
				inln += "\t" + currentFDR;
			}
			
			if (currentFDR < threshold) {
				significant++;
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
				System.out.print(ctr + " lines parsed... " + significant + "significant\r");
			}
			inln = tf.readLine();
		}
		
		System.out.println();
		System.out.println();
		System.out.print("Done. w" + ctr + " lines parsed in total. " + significant + "significant\r");
		
		tfosig.close();
		if (!onlywritesignificant) {
			tfo.close();
		}
		tf.close();
	}
	
}
