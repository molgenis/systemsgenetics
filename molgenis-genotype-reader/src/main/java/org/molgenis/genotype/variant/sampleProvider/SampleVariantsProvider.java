package org.molgenis.genotype.variant.sampleProvider;

import java.util.List;

import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 * Loads the sample variants for a variant. Is used to enable lazy loading of
 * the sample variants
 * 
 * @author erwin
 * 
 */
public interface SampleVariantsProvider
{
	/**
	 * Get sample variants. Do not call directly always access via the
	 * variant.getSampleVariants
	 * 
	 * @param variant
	 * @return
	 */
	List<Alleles> getSampleVariants(GeneticVariant variant);

	/**
	 * Returns for each sample if it phased or not
	 * 
	 * @return
	 */
	List<Boolean> getSamplePhasing(GeneticVariant variant);

	int cacheSize();

	int getSampleVariantProviderUniqueId();

	/**
	 * Get sample called dosage {0,1,2} -1 denotes missing
	 * 
	 * @return
	 */
	byte[] getSampleCalledDosage(GeneticVariant variant);

	/**
	 * Get sample dosage in range of 0 to 2. -1 denotes missing
	 * 
	 * @return
	 */
	float[] getSampleDosage(GeneticVariant variant);

}
