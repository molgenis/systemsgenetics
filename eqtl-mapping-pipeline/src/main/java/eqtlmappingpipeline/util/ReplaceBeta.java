package eqtlmappingpipeline.util;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashMap;

public class ReplaceBeta {
	
	public void run(String in, String out, String mafqcfile) throws IOException {
		TextFile tf = new TextFile(in, TextFile.R);
		TextFile tfo = new TextFile(out, TextFile.W);
		String[] header = tf.readLineElems(TextFile.tab);
		
		int metazcol = -1;
		int ncol = -1;
		int metabetacol = -1;
		int snpcol = -1;
		
		for (int i = 0; i < header.length; i++) {
			
			if (header[i].contains("")) {
				metazcol = i;
			} else if (header[i].contains("")) {
				metabetacol = i;
			} else if (header[i].contains("")) {
				snpcol = i;
			} else if (header[i].contains("")) {
				ncol = i;
			}
			
			
		}
		
		if (metazcol < 0 || metabetacol < 0 || ncol < 0 || snpcol < 0) {
			System.out.println("Error: could not find one of the columns. Found the following:");
			System.out.println("MetaZ: " + metazcol);
			System.out.println("SNP: " + snpcol);
			System.out.println("MetaBeta: " + metabetacol);
			System.out.println("Ncol: " + ncol);
			System.exit(-1);
		}
		
		// load maf per snp
		HashMap<String, Double> mafpersnp = new HashMap<String, Double>();
		TextFile tf1 = new TextFile(mafqcfile, TextFile.R);
		tf1.readLine();
		String[] elems = tf1.readLineElems(TextFile.tab);
		while (elems != null) {
			String snp = elems[0];
			double maf = Double.parseDouble(elems[elems.length - 1]);
			mafpersnp.put(snp, maf);
			elems = tf1.readLineElems(TextFile.tab);
		}
		tf1.close();
		
		elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			
			String snp = elems[snpcol];
			
			double maf = mafpersnp.get(snp);
			String nstr = elems[ncol];
			int n = 0;
			String[] nstrelems = nstr.split(";");
			for (String s : nstrelems) {
				if (!s.equals("-")) {
					n += Integer.parseInt(s);
				}
			}
			double z = Double.parseDouble(elems[metazcol]);
			
			double[] betaandse = ZScores.zToBeta(z, maf, n);
			elems[metabetacol] = betaandse[0] + " (" + betaandse[1] + ")";
			
			tfo.writeln(Strings.concat(elems, Strings.tab));
			elems = tf.readLineElems(TextFile.tab);
		}
		
		tf.close();
		tfo.close();
	}
	
}
