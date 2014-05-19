package org.molgenis.genotype.variant.sampleProvider;

import java.util.List;

import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.util.FixedSizeIterable;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.GenotypeRecord;

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
	 * Get information of each sample variant for a variant
	 * 
	 * @param variant
	 * @return
	 */
	FixedSizeIterable<GenotypeRecord> getSampleGenotypeRecords(GeneticVariant variant);

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
	
	/**
	 * Get sample genotype probabilities. Make sure to ask genotype data if this is save
	 * 
	 * [sample][AA,AB,BB]
	 * 
	 * @param variant
	 * @return 
	 */
	float[][] getSampleProbilities(GeneticVariant variant);

}
