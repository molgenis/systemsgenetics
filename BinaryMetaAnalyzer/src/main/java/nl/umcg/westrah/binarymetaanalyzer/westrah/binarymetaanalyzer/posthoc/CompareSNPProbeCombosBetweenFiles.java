package nl.umcg.westrah.binarymetaanalyzer.westrah.binarymetaanalyzer.posthoc;

import umcg.genetica.io.text.TextFile;

import javax.xml.soap.Text;
import java.io.IOException;
import java.util.HashSet;

public class CompareSNPProbeCombosBetweenFiles {
	
	public void run(String f1, String f2, String fout) throws IOException {
		TextFile tf = new TextFile(f1, TextFile.R, 10 * 1048576);
		System.out.println(f1);
		tf.readLine();
		HashSet<String> set1 = new HashSet<String>(140000000);
		HashSet<String> set2 = new HashSet<String>(140000000);
		String[] elems = tf.readLineElems(TextFile.tab);
		int ctr = 0;
		while (elems != null) {
			
			String snp = elems[1];
			String gene = elems[2];
			
			String meh = snp + "__" + gene;
			set1.add(meh);
			ctr++;
			if (ctr % 10000 == 0) {
				System.out.print(ctr + " read.\r");
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		System.out.println();
		
		TextFile tf2 = new TextFile(f2, TextFile.R, 10 * 1048576);
		System.out.println(f2);
		tf2.readLine();
		elems = tf2.readLineElems(TextFile.tab);
		ctr = 0;
		while (elems != null) {
			String snp = elems[1];
			String gene = elems[2];
			
			String meh = snp + "__" + gene;
			set2.add(meh);
			ctr++;
			if (ctr % 10000 == 0) {
				System.out.print(ctr + " read.\r");
			}
			elems = tf2.readLineElems(TextFile.tab);
		}
		System.out.println();
		
		tf2.close();
		
		
		int shared = 0;
		int notfound = 0;
		TextFile outf = new TextFile(fout, TextFile.W);
		for (String str : set1) {
			if (!set2.contains(str)) {
				String[] elmes = str.split("__");
				outf.writeln(elmes[0] + "\t" + elmes[1]);
				notfound++;
			} else {
				shared++;
			}
		}
		outf.close();
		System.out.println(notfound + " not found");
		System.out.println(shared + " shared");
		
		
	}
}
