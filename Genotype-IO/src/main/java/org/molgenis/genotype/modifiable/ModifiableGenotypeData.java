package org.molgenis.genotype.modifiable;

import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.id.GeneticVariantId;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;

/**
 * GenotypeData of which some variants properties can be updated and from which
 * variants can be excluded
 * 
 * @author Patrick Deelen
 * 
 */
public interface ModifiableGenotypeData extends RandomAccessGenotypeData
{

	/**
	 * Get the updated ID for a genetic variant. Return null if not updated.
	 * 
	 * @param geneticVariant
	 *            the potentially modified variant
	 * @return the updated ID or null if not modified
	 */
	GeneticVariantId getUpdatedId(ModifiableGeneticVariant geneticVariant);

	/**
	 * Get the updated reference allele for a genetic variant. Return null if
	 * not updated.
	 * 
	 * @param geneticVariant
	 *            the potentially modified variant
	 * @return the updated reference allele or null if not modified
	 */
	Allele getUpdatedRef(ModifiableGeneticVariant geneticVariant);

	/**
	 * Get the updated sample variant provider
	 * 
	 * @param geneticVariant
	 *            the potentially modified variant
	 * @return the updated sample provider or null if not modified
	 */
	SampleVariantsProvider getUpdatedSampleVariantProvider(ModifiableGeneticVariant geneticVariant);

	/**
	 * Overwrite the original genetic variant ID.
	 * 
	 * @param geneticVariant
	 *            the variant to modify
	 * @param geneticVariantId
	 *            the new genetic variant ID or null if not modified
	 */
	void updateVariantId(ModifiableGeneticVariant geneticVariant, GeneticVariantId newGeneticVariantId);

	/**
	 * Updates the original genetic variant ID to make the new key primary. Old
	 * keys will be made alternative IDs.
	 * 
	 * @param geneticVariant
	 *            the variant to modify
	 * @param newPrimaryId
	 *            the new primary ID or null if not modified
	 */
	void updateVariantPrimaryId(ModifiableGeneticVariant geneticVariant, String newPrimaryId);

	/**
	 * Updates the the sample variant provider and the ref allele. The sample
	 * variant provider will be the swapping sample variant provider and the
	 * reference allele will be made complement
	 * 
	 * @param geneticVariant
	 *            the variant to modify
	 */
	void swapGeneticVariant(ModifiableGeneticVariant geneticVariant);

	/**
	 * Updates the reference allele. Also changes the order of the alleles to
	 * make reference allele first
	 * 
	 * @param geneticVariant
	 *            the variant to modify
	 * @param newRefAllele
	 * @throws GenotypeModificationException
	 *             if reference allele is not one of the variant alleles
	 */
	void updateRefAllele(ModifiableGeneticVariant geneticVariant, Allele newRefAllele)
			throws GenotypeModificationException;

	/**
	 * Get the updated alleles
	 * 
	 * @param geneticVariant
	 *            the potentially modified variant
	 * @return the updated alleles or null if not modified
	 */
	Alleles getUpdatedAlleles(ModifiableGeneticVariant geneticVariant);

	/**
	 * Get all modifiable variants from a sequence
	 * 
	 * @param seqName
	 * @return
	 */
	Iterable<ModifiableGeneticVariant> getModifiableSequenceGeneticVariants(String seqName);

	/**
	 * Get the modifiable variants at the specified position
	 * 
	 * @param seqName
	 * @param startPos
	 * @return all variants found at startPos, can be empty if none found
	 */
	Iterable<ModifiableGeneticVariant> getModifiableVariantsByPos(String seqName, int startPos);

	/**
	 * Get the modifiable SNP variant at the specified position. Only one SNP
	 * possible per position.
	 * 
	 * @param seqName
	 * @param startPos
	 * @return The SNP found at this startPos, will be null if not present
	 */
	ModifiableGeneticVariant getModifiableSnpVariantByPos(String seqName, int startPos);

	/**
	 * Get all modifiable variants
	 * 
	 * @return
	 */
	Iterable<ModifiableGeneticVariant> getModifiableGeneticVariants();

	/**
	 * Exclude a variant from future queries on this dataset
	 * 
	 * @param geneticVariant
	 */
	void excludeVariant(ModifiableGeneticVariant geneticVariant);

	/**
	 * Get the number of excluded variants
	 * 
	 * @return excluded variant count
	 */
	int getExcludedVariantCount();
	
	public boolean isSwapped(GeneticVariant geneticVariant);

}
