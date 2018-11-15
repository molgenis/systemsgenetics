package nl.umcg.westrah.binarymetaanalyzer.westrah.binarymetaanalyzer.posthoc;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;

import java.io.IOException;
import java.util.Iterator;

public class BHFDR {
	
	public static void main(String[] args) {
		String in = "D:\\biogen\\trans\\eQTLs.txt.gz";
		String out = "D:\\biogen\\trans\\eQTLs-fdr.txt.gz";
		double p = 0.05;
		BHFDR f = new BHFDR();
		try {
			f.run(in, p, out);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	public void run(String in, double significance, String out) throws IOException {
		System.out.println("Benjamini Hochberg correction");
		System.out.println("in: " + in);
		System.out.println("out: " + out);
		System.out.println("Threshold: " + significance);
		TextFile tf = new TextFile(in, TextFile.R);
		int nrlines = tf.countLines();
		tf.close();
		
		int nrtests = nrlines - 1;
		QTLTextFile tf2 = new QTLTextFile(in, TextFile.R);
		TextFile outf = new TextFile(out, TextFile.W);
		TextFile outfsig = new TextFile(out + "-FDR" + significance + ".txt.gz", TextFile.W);
		String header = tf2.readLine();
		outf.writeln(header);
		outfsig.writeln(header);
		Iterator<EQTL> it = tf2.getEQtlIterator();
		
		int ctr = 1;
		double prevp = -1;
		int nrsig = 0;
		
		while (it.hasNext()) {
			EQTL e = it.next();
			if (prevp == -1) {
				prevp = e.getPvalue();
			} else {
				if (e.getPvalue() < prevp) {
					System.out.println("QTL File not properly sorted! Faulty line: " + ctr + "\tSNP: " + e.getRsName() + "\tGene: " + e.getProbe());
					System.exit(-1);
				}
				prevp = e.getPvalue();
			}
			double qvalue = e.getPvalue() * (nrtests / ctr);
//			double qvalue = ((ctr + 1d) / nrlines) * e.getPvalue();
//			System.out.println(ctr + "\t" + qvalue);
//			if (ctr > 10) {
//				System.exit(-1);
//			}
			e.setFDR(qvalue);
			outf.writeln(e.toString());
			if (qvalue < significance) {
				outfsig.writeln(e.toString());
				nrsig++;
			}
			ctr++;
		}
		tf2.close();
		outf.close();
		outfsig.close();
		System.out.println(nrsig + " out of " + nrtests + " are signifciant.");
	}
	
}
