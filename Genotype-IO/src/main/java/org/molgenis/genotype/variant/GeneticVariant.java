package org.molgenis.genotype.variant;

import java.util.List;
import java.util.Map;

import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.util.FixedSizeIterable;
import org.molgenis.genotype.util.Ld;
import org.molgenis.genotype.util.LdCalculatorException;
import org.molgenis.genotype.variant.id.GeneticVariantId;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;

public interface GeneticVariant extends Comparable<GeneticVariant>
{
	/**
	 * Get the meta data for this variant
	 * 
	 * @return
	 */
	public GeneticVariantMeta getVariantMeta();
	
	/**
	 * Get the primary variant ID
	 * 
	 * @return String
	 */
	public String getPrimaryVariantId();

	/**
	 * Gets all the other id's (names) besides the primaryVariantId by which
	 * this variant is known.
	 * 
	 * @return List of String
	 */
	public List<String> getAlternativeVariantIds();

	/**
	 * Get all IDs for this variant
	 * 
	 * @return List of String
	 */
	public List<String> getAllIds();

	/**
	 * Get the variant ID object for this variant
	 * 
	 * @return
	 */
	public GeneticVariantId getVariantId();

	/**
	 * Gets the starting position on the sequence
	 * 
	 * @return int
	 */
	public int getStartPos();

	/**
	 * Get the Sequence this variant is located on
	 * 
	 * @return the Sequence
	 */
	public String getSequenceName();

	/**
	 * Get all possible alleles (including the reference) The first value is the
	 * reference value if it is set
	 * 
	 * @return
	 */
	public Alleles getVariantAlleles();

	/**
	 * Get the total allele count
	 * 
	 * @return
	 */
	public int getAlleleCount();

	/**
	 * Gets the reference allele
	 * 
	 * @return String
	 */
	public Allele getRefAllele();
	
	/**
	 * Gets the alternative alleles
	 * 
	 */
	public Alleles getAlternativeAlleles();

	/**
	 * Returns list sample variants.
	 */
	public List<Alleles> getSampleVariants();

	/**
	 * Get the annotations for this variant. The key is the annotationId, the
	 * value is of the type defined by the Annotation (see
	 * GenotypeData.getVariantAnnotations()
	 * 
	 * @return
	 */
	public Map<String, ?> getAnnotationValues();

	/**
	 * Get the frequency of the minor allele
	 * 
	 * @return the minor allele frequency
	 */
	public double getMinorAlleleFrequency();

	/**
	 * Get the minor allele
	 * 
	 * @return the minor allele
	 */
	public Allele getMinorAllele();

	/**
	 * Is this variant a SNP
	 * 
	 * @return true if SNP
	 */
	public boolean isSnp();

	/**
	 * Is this variant an AT or GC SNP.
	 * 
	 * @return true if SNP and AT or GC
	 */
	public boolean isAtOrGcSnp();

	/**
	 * Calculate LD to an other variant
	 * 
	 * @param other
	 * @return
	 * @throws LdCalculatorException
	 */
	public Ld calculateLd(GeneticVariant other) throws LdCalculatorException;

	/**
	 * Test if biallelic
	 * 
	 * @return
	 */
	public boolean isBiallelic();

	/**
	 * Get the dosage values. -1 for unknown
	 * 
	 * @return dosage values
	 */
	public float[] getSampleDosages();

	/**
	 * Dosage values of 0, 1 or 2. Count of reference allele. -1 for missing
	 * data. If ref is null then count first of variant alleles.
	 * 
	 * @return
	 */
	public byte[] getSampleCalledDosages();

	/**
	 * Gets boolean values for every sample indicating if, for this variant,
	 * phased data is present in the sample.
	 *
	 * @return a list of booleans indicating presence of phasing data for every sample.
	 */
	public List<Boolean> getSamplePhasing();

	/**
	 * Method returning if phased data is present or not.
	 *
	 * @return true if every sample contains phased data for this
	 * variant, otherwise false.
	 */
	public boolean hasPhasedGenotypes();

	/**
	 * Method returning if phased data is present or not.
	 *
	 * @return true if probabilities are available for every sample for this variant,
	 * otherwise false.
	 */
	public boolean hasPhasedProbabilities();

	/**
	 * Get the sample variant provider used by this variant
	 * 
	 * @return
	 */
	public SampleVariantsProvider getSampleVariantsProvider();

	public boolean isMapped();

	public double getCallRate();

	/**
	 * Warning currently only HWE calculation for bi-allelic variants
	 * Other variants will have HWE of NaN
	 * 
	 * @return 
	 */
	public double getHwePvalue();
	
	/**
	 * Calculate the MACH r2 measure
	 * 
	 * For formula see: doi:10.1038/nrg2796 S3
	 * 
	 * @return 
	 */
	public double getMachR2();
	
	/**
	 * [sample][AA,AB,BB]
	 * 
	 * @return 
	 */
	public float[][] getSampleGenotypeProbilities();

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
	 * @return An array of probabilities for every sample and possible genotype for a sample and the given variant
	 */
	public double[][] getSampleGenotypeProbabilitiesComplex();

	/**
	 * [sample][haplotype][A, B, C, ...]
	 *
	 * Get sample haplotype probabilities. For every allele per haplotype a probability is stored for a particular sample.
	 * Make sure to ask whether phased data is available for this variant.
	 *
	 * @return An array of probabilities per haplotype, per sample.
	 */
	public double[][][] getSampleGenotypeProbabilitiesPhased();

	/**
	 * Get the records that are available for this variants from the different samples
	 * 
	 * @return 
	 */
	public FixedSizeIterable<GenotypeRecord> getSampleGenotypeRecords();

}
