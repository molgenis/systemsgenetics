package umcg.genetica.io.binInteraction.gene;

/**
 *
 * @author Patrick Deelen
 */
public interface BinaryInteractionGene {

	String getName();

	String getChr();

	int getStart();

	int getEnd();
	
	int[] getVariantPointers();
	
	int getVariantCount();
	
	int getIndexOfVariantPointer(int variantIndex);
	
}
