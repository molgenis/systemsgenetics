package umcg.genetica.features;

import umcg.genetica.enums.Chromosome;
import umcg.genetica.enums.Strand;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by hwestra on 9/13/16.
 */
public class GFFFile {

	public ArrayList<Feature> readGFF(String gff, boolean onlySuggestedLoci) throws IOException {
		TextFile tf = new TextFile(gff, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		ArrayList<Feature> output = new ArrayList<Feature>();
		while (elems != null) {
			if (elems[0].startsWith("#")) {

			} else {

				Chromosome chr = Chromosome.parseChr(elems[0]);
				Integer start = Integer.parseInt(elems[3]);
				Integer stop = Integer.parseInt(elems[4]);
				String info = elems[8];
				String[] infoElems = info.split(";");
				boolean isCandidate = false;
				Feature f = new Feature();
				f.setChromosome(chr);
				f.setStrand(Strand.NEG);
				f.setStart(start);
				f.setStop(stop);
				for (int i = 0; i < infoElems.length; i++) {
					if (infoElems[i].startsWith("is_candidate")) {
						String[] data = infoElems[i].split("=");
						if (data[1].equals("1")) {
							isCandidate = true;
						}
					} else if (infoElems[i].startsWith("Name")) {
						String[] data = infoElems[i].split("=");
						f.setName(data[1]);
					}
				}


				if (onlySuggestedLoci) {
					if (isCandidate) {
						if (f.getName().equals("ACOXL")) {
							System.out.println("found it");
						}
						output.add(f);
					}
				} else {
					output.add(f);
				}


			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		return output;
	}
}
