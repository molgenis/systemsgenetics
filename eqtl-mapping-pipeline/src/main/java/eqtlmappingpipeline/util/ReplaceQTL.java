package eqtlmappingpipeline.util;

import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashMap;

public class ReplaceQTL {
	
	
	public static void main(String[] args) {
		String eqtl = "D:\\Sync\\SyncThing\\Postdoc2\\2017-11-eQTLMeta\\data\\2018-04-10-Replication\\2018-04-03-eQTLsFDR-Significant-0.05-WithoutCrossMappingEffects-PrunedFDR-Significant-0.05.txt.gz";
		String eqtm = "D:\\Sync\\SyncThing\\Postdoc2\\2017-11-eQTLMeta\\data\\2018-04-10-Replication\\2018-04-10-cis-eqtm-1mb.txt.gz";
		String out = "D:\\Sync\\SyncThing\\Postdoc2\\2017-11-eQTLMeta\\data\\2018-04-10-Replication\\2018-04-11-eQTLsFDR-Significant-0.05-WithoutCrossMappingEffects-PrunedFDR-Significant-0.05-CpG.txt.gz";
		
		ReplaceQTL q = new ReplaceQTL();

//		try {
//			q.convertUsingLinkFile(eqtm, eqtl, out, false, false);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
		
		String meqtl = "D:\\Sync\\SyncThing\\Postdoc2\\2017-11-eQTLMeta\\data\\2018-04-10-Replication\\2017-04-10-trans-meqtl-transrepl.txt.gz";
		out = "D:\\Sync\\SyncThing\\Postdoc2\\2017-11-eQTLMeta\\data\\2018-04-10-Replication\\2017-04-10-trans-meqtl-transrepl-snpgenepairs.txt.gz";
		try {
			q.convertUsingLinkFile(eqtm, meqtl, out, false, true);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		meqtl = "D:\\Sync\\SyncThing\\Postdoc2\\2017-11-eQTLMeta\\data\\2018-04-10-Replication\\2017-04-10-trans-meqtl-transrepl-FDR0.05.txt.gz";
		out = "D:\\Sync\\SyncThing\\Postdoc2\\2017-11-eQTLMeta\\data\\2018-04-10-Replication\\2017-04-10-trans-meqtl-transrepl-snpgenepairs-FDR0.05.txt.gz";
		try {
			q.convertUsingLinkFile(eqtm, meqtl, out, false, true);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	
	public void convertUsingLinkFile(String linkfile, String input, String out, boolean alsowritenonmappingmeqtls, boolean aToB) throws IOException {
		HashMap<String, Pair<String, Integer>> replacements = new HashMap<>();
		
		
		TextFile tf = new TextFile(linkfile, TextFile.R);
		System.out.println("Reading: " + linkfile);
		String[] header = tf.readLineElems(TextFile.tab);
		boolean eqtlfile = false;
		if (header.length > 20) {
			System.out.println("Header suggests eQTL file for linkfile");
			eqtlfile = true;
		} else if (header.length == 3) {
			System.out.println("Assuming linkfile is not an eQTL file");
		} else {
			System.out.println("Please input linkfile with > 3 columns");
			tf.close();
			System.exit(-1);
		}
		
		int lnctr = 0;
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String a = null;
			String b = null;
			Integer d = 1;
			
			if (!eqtlfile) {
				a = new String(elems[0].getBytes("UTF-8")).intern();
				b = new String(elems[1].getBytes("UTF-8")).intern();
				String z = elems[2];
				if (Double.parseDouble(z) < 0) {
					d = -1;
				}
			} else {
				EQTL e = EQTL.fromString(elems, "-", TextFile.semicolon);
				
				a = e.getRsName();
				b = e.getProbe();
				
				// get the alleles, flip if required.
				double z = e.getZscore();
				
				if (!e.getAlleleAssessed().equals("C")) {
					z *= -1;
				}
				
				if (z < 0) {
					d = -1;
				}
			}
			
			if (aToB) {
				Pair<String, Integer> q = new Pair<>(b, d);
				if (!replacements.containsKey(a)) {
					replacements.put(a, q);
				}
				
			} else {
				Pair<String, Integer> q = new Pair<>(a, d);
				if (!replacements.containsKey(b)) {
					replacements.put(b, q);
				}
			}
			
			elems = tf.readLineElems(TextFile.tab);
			lnctr++;
		}
		tf.close();
		
		System.out.println(lnctr + " lines in " + linkfile);
		System.out.println(replacements.size() + " replacements loaded");
		
		TextFile tf2 = new TextFile(input, TextFile.R);
		TextFile tf3 = new TextFile(out, TextFile.W);
		
		System.out.println("Replacing QTLs in file: " + input);
		
		tf3.writeln(tf2.readLine());
		
		int replaced = 0;
		elems = tf2.readLineElems(TextFile.tab);
		while (elems != null) {
			if (elems.length > 4) {
				String snp = elems[1];
				String probe = elems[4];
				
				
				Pair<String, Integer> d = replacements.get(probe);
				if (d != null) {
					elems[4] = d.getLeft();
					Double z = Double.parseDouble(elems[10]);
					z *= d.getRight();
					elems[10] = "" + z;
					tf3.writeln(Strings.concat(elems, Strings.tab));
					replaced++;
				} else if (alsowritenonmappingmeqtls) {
					tf3.writeln(Strings.concat(elems, Strings.tab));
				}
				
			}
			elems = tf2.readLineElems(TextFile.tab);
		}
		
		tf3.close();
		tf2.close();
		
		System.out.println(replaced + " lines replaced.");
		
	}
	
	
}
