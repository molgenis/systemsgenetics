package eqtlmappingpipeline.util;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.ProbeAnnotation;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

public class BatchSplitter {
	
	public void splitPhenotypePerChr(String matrix, String probeAnnotation, String outLoc) throws IOException {
		
		ProbeAnnotation pb = new ProbeAnnotation();
		System.out.println("Loading: " + probeAnnotation);
		pb.load(probeAnnotation);
		
		HashSet<Integer> uniqueChrs = new HashSet<>();
		String[] probes = pb.getProbes();
		for (int i = 0; i < probes.length; i++) {
			if (pb.getChr()[i] > 0) {
				uniqueChrs.add((int) pb.getChr()[i]);
			}
		}
		
		
		TextFile tf = new TextFile(matrix, TextFile.R);
		String header = tf.readLine(); // header
		TextFile[] out = new TextFile[uniqueChrs.size()];
		HashMap<Integer, Integer> chrToFile = new HashMap<Integer, Integer>();
		int ctr = 0;
		for (Integer key : uniqueChrs) {
			chrToFile.put(key, ctr);
			String filename = outLoc + "-chr" + key + ".txt.gz";
			System.out.println("Creating file: " + filename);
			out[ctr] = new TextFile(filename, TextFile.W);
			
			out[ctr].writeln(header);
			ctr++;
		}
		
		String ln = tf.readLine();
		HashMap<String, Integer> probeIds = pb.getProbeToProbeId();
		short[] chrs = pb.getChr();
		int parsed = 0;
		int written = 0;
		while (ln != null) {
			String[] elems = null;
			if (ln.length() > 200) {
				elems = ln.substring(0, 200).split("\t");
			} else {
				elems = ln.split("\t");
			}
			
			String probe = elems[0];
			Integer probeid = probeIds.get(probe);
			if (probeid != null) {
				int chr = (int) chrs[probeid];
				if (chr > 0) {
					Integer fileId = chrToFile.get(chr);
					out[fileId].writeln(ln);
					written++;
				}
			}
			parsed++;
			
			if (parsed % 1000 == 0) {
				System.out.print("\rLines parsed: " + parsed + "\twritten: " + written);
			}
			ln = tf.readLine();
		}
		
		for (TextFile t : out) {
			if (t != null) {
				t.close();
			}
		}
		System.out.println("Done");
		
	}
	
}
