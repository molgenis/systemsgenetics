package nl.umcg.westrah.binarymetaanalyzer.westrah.binarymetaanalyzer.posthoc;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class SplitQTLFile {
	
	public static void main(String[] args) {
	
	}
	
	public void split(String efile, String groupdefinition, String output) throws IOException {
		
		
		TextFile tf = new TextFile(groupdefinition, TextFile.R);
		String[] elms = tf.readLineElems(TextFile.tab);
		
		HashMap<String, Integer> groupIndex = new HashMap<String, Integer>();
		ArrayList<String> groupList = new ArrayList<String>();
		ArrayList<HashSet<String>> grpDsSets = new ArrayList<>();
		while (elms != null) {
			if (elms.length > 1) {
				String ds = elms[0];
				String grp = elms[1];
				
				Integer grpid = groupIndex.get(grp);
				if (grpid == null) {
					groupIndex.put(grp, groupList.size());
					groupList.add(grp);
					grpDsSets.add(new HashSet<>());
				}
				
				HashSet<String> grpSet = grpDsSets.get(grpid);
				grpSet.add(ds);
				
			}
			elms = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		TextFile[] outputs = new TextFile[groupList.size()];
		TextFile input = new TextFile(efile, TextFile.R);
		String header = input.readLine();
		for (int t = 0; t < outputs.length; t++) {
			String grp = groupList.get(t);
			outputs[t] = new TextFile(output + "-" + grp + ".txt.gz", TextFile.W);
			outputs[t].writeln(header);
		}
		
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String[] datasets = elems[10].split(";");
			String[] datasetsZ = elems[10].split(";");
			String[] datasetsN = elems[11].split(";");
			
			for (int g = 0; g < groupList.size(); g++) {
				ArrayList<Double> z = new ArrayList<>();
				ArrayList<Integer> n = new ArrayList<>();
				ArrayList<String> ds = new ArrayList<>();
				
				HashSet<String> grpSet = grpDsSets.get(g);
				for (int d = 0; d < datasets.length; d++) {
					if (grpSet.contains(datasets[d])) {
						z.add(Double.parseDouble(datasetsZ[d]));
						n.add(Integer.parseInt(datasetsN[d]));
						ds.add(datasets[d]);
					}
				}
				
				double[] dz = Primitives.toPrimitiveArr(z);
				int[] in = Primitives.toPrimitiveArr(n);
				double metaz = ZScores.getWeightedZ(dz, in);
				
				String[] outputE = new String[elems.length];
				System.arraycopy(elems, 0, outputE, 0, elems.length);
				outputE[10] = "" + metaz;
				outputE[11] = Strings.concat(dz, Strings.semicolon);
				outputE[12] = Strings.concat(in, Strings.semicolon);
				outputE[13] = Strings.concat(ds, Strings.semicolon);
				outputs[g].writeln(Strings.concat(outputE, Strings.tab));
			}
			
			elems = tf.readLineElems(TextFile.tab);
		}
		
		
		for (int t = 0; t < outputs.length; t++) {
			outputs[t].close();
		}
	}
	
}
