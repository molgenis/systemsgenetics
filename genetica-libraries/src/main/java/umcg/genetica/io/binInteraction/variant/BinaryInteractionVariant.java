package umcg.genetica.io.binInteraction.variant;

import org.molgenis.genotype.Allele;

/**
 *
 * @author Patrick Deelen
 */
public interface BinaryInteractionVariant {

	Allele getAltAllele();

	String getChr();

	int getGeneCount();

	int[] getGenePointers();

	String getName();

	int getPos();

	Allele getRefAllele();
	
	int getIndexOfGenePointer(int geneIndex);

}
