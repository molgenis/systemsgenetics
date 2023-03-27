package mbqtl.datastructures;

import umcg.genetica.enums.Chromosome;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class SNPAnnotation {
	HashMap<String, Integer> variantToId = new HashMap<>();
	ArrayList<String> variants = new ArrayList<>();
	ArrayList<Chromosome> chromosomes = new ArrayList<>();
	ArrayList<Integer> positions = new ArrayList<>();

	public SNPAnnotation(String snpannotation) throws IOException {

		load(snpannotation);
	}

	private void load(String snpannotation) throws IOException {
		System.out.println("Loading variant annotation from: " + snpannotation);
		TextFile tf = new TextFile(snpannotation, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		int ctr = 0;
		while (elems != null) {
			if (elems.length > 2) {
				try {
					String variant = elems[0];
					Chromosome chr = Chromosome.parseChr(elems[1]);
					int pos = Integer.parseInt(elems[2]);

					if (variantToId.containsKey(variant)) {
						int tmpid = variantToId.get(variant);
						Chromosome tmpchr = chromosomes.get(tmpid);
						Integer tmppos = positions.get(tmpid);
						if (!tmpchr.equals(chr) || !tmppos.equals(pos)) {
							System.out.println("Variant " + variant + " has multiple annotations. Is this correct? Only storing first instance.");
						}
					} else {
						variantToId.put(variant, ctr);
						variants.add(variant);
						chromosomes.add(chr);
						positions.add(pos);
						ctr++;
					}
				} catch (NumberFormatException e) {

				}
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		System.out.println("Loaded annotation for " + variants.size() + " variants.");
	}

	public Integer getId(String variantName) {
		return variantToId.get(variantName);
	}

	public int getPos(Integer id) {
		return positions.get(id);
	}
}
