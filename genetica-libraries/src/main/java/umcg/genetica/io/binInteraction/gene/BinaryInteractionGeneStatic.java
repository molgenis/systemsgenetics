package umcg.genetica.io.binInteraction.gene;

/**
 *
 * @author Patrick Deelen
 */
public class BinaryInteractionGeneStatic extends BinaryInteractionGeneAbstract {
	
	private final int[] variants;

	public BinaryInteractionGeneStatic(String name, String chr, int start, int end, int[] variants) {
		super(name, chr, start, end);
		this.variants = variants;
	}

	@Override
	public int[] getVariantPointers() {
		return variants;
	}

	@Override
	public int getVariantCount() {
		return variants.length;
	}

}
