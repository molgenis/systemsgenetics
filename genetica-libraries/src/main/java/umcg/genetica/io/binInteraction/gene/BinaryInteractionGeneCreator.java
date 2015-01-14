package umcg.genetica.io.binInteraction.gene;

import gnu.trove.list.array.TIntArrayList;

/**
 *
 * @author Patrick Deelen
 */
public class BinaryInteractionGeneCreator extends BinaryInteractionGeneAbstract {
	
	private TIntArrayList variantPointers;

	public BinaryInteractionGeneCreator(String name, String chr, int start, int end) {
		super(name, chr, start, end);
		variantPointers = new TIntArrayList();
	}

	@Override
	public int[] getVariantPointers() {
		return variantPointers.toArray();
	}

	@Override
	public int getVariantCount() {
		return variantPointers.size();
	}
	
	public void addVariant(int variantPointer){
		variantPointers.add(variantPointer);
	}

}
