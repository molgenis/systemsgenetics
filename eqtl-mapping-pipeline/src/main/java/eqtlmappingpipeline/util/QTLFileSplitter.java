package eqtlmappingpipeline.util;


import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashMap;

public class QTLFileSplitter {
	
	
	public void split(String in, String out) throws IOException {
		
		TextFile tf = new TextFile(in, TextFile.R);
		String header = tf.readLine();
		String ln = tf.readLine();
		
		HashMap<String, TextFile> chrmap = new HashMap<String, TextFile>();
		int nrlines = 0;
		while (ln != null) {
			String[] elems = Strings.subsplit(ln, Strings.tab, 2, 3);
			TextFile chrout = chrmap.get(elems[0]);
			if (chrout == null) {
				chrout = new TextFile(out + "-" + elems[0] + ".txt.gz", TextFile.W);
				System.out.println("Opening outputfile: " + chrout.getFileName());
				chrout.writeln(header);
				chrmap.put(elems[0], chrout);
			}
			chrout.writeln(ln);
			nrlines++;
			if (nrlines % 10000 == 0) {
				System.out.print("\r" + nrlines + "\tprocessed.");
			}
			ln = tf.readLine();
		}
		tf.close();
		System.out.println("");
		
		for (String key : chrmap.keySet()) {
			chrmap.get(key).close();
		}
		System.out.println("done.");
	}
	
}
