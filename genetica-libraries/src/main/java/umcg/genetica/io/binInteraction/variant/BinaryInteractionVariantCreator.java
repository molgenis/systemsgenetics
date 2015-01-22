package umcg.genetica.io.binInteraction.variant;

import gnu.trove.list.array.TIntArrayList;
import org.molgenis.genotype.Allele;

/**
 *
 * @author Patrick Deelen
 */
public class BinaryInteractionVariantCreator extends BinaryInteractionVariantAbstract {

	private TIntArrayList genePointers;

	public BinaryInteractionVariantCreator(String name, String chr, int pos, Allele refAllele, Allele altAllele) {
		super(name, chr, pos, refAllele, altAllele);
		genePointers = new TIntArrayList();
	}

	@Override
	public int getGeneCount() {
		return genePointers.size();
	}

	@Override
	public int[] getGenePointers() {
		return genePointers.toArray();
	}

	public void addGene(int genePointer) {
		genePointers.add(genePointer);
	}
	
	public void sortGenePointers(){
		genePointers.sort();
	}

	@Override
	public int getIndexOfGenePointer(int geneIndex) {
		return genePointers.binarySearch(geneIndex);
	}
	
}
