package umcg.genetica.io.binInteraction.variant;

import java.util.Arrays;
import org.molgenis.genotype.Allele;

/**
 *
 * @author Patrick Deelen
 */
public class BinaryInteractionVariantStatic extends BinaryInteractionVariantAbstract {
	
	private final int[] genePointers;

	public BinaryInteractionVariantStatic(String name, String chr, int pos, Allele refAllele, Allele altAllele, int[] genePointers) {
		super(name, chr, pos, refAllele, altAllele);
		this.genePointers = genePointers;
	}

	@Override
	public int getGeneCount() {
		return genePointers.length;
	}

	@Override
	public int[] getGenePointers() {
		return genePointers;
	}

	@Override
	public int getIndexOfGenePointer(int geneIndex) {
		return Arrays.binarySearch(genePointers, geneIndex);
	}
	

}
