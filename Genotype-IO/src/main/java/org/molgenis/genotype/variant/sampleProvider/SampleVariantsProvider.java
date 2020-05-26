package org.molgenis.genotype.variant.sampleProvider;

import java.util.List;

import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.util.FixedSizeIterable;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.GenotypeRecord;
import org.molgenis.genotype.variant.ReadOnlyGeneticVariantBgen;

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

	/**
	 * [sample][AAA,AAB,ABB,...,ACC,BCC,CCC]
	 *
	 * Following <a href="https://www.well.ox.ac.uk/~gav/bgen_format/spec/latest.html">BGEN specification</a>,
	 * the probabilities per sample are stored for every possible genotype given
	 * the ploidy of the sample and the possible alleles. Probabilities are stored in colexicographic order of the
	 * K-vectors of nonnegative integers (X1, X2, ..., Xk) representing the count of the i-th allele in the genotype.
	 *
	 * In contrast to the BGEN specification, all probabilities are stored so the sum of these is one.
	 *
	 * @param variant The variant to request probabilities for.
	 * @return An array of probabilities for every sample and possible genotype for a sample and the given variant
	 */
	public double[][] getSampleProbabilitiesComplex(GeneticVariant variant);

	/**
	 * [sample][haplotype][A, B, C, ...]
	 *
	 * Get sample haplotype probabilities. For every allele per haplotype a probability is stored for a particular sample.
	 * Make sure to ask whether phased data is available for the variant.
	 *
	 * @param variant The variant to request probabilities for.
	 * @return An array of probabilities per haplotype, per sample.
	 */
	public double[][][] getSampleProbabilitiesPhased(GeneticVariant variant);
}
