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

	public List<Boolean> getSamplePhasing();

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
	 * Get the records that are available for this variants from the different samples
	 * 
	 * @return 
	 */
	public FixedSizeIterable<GenotypeRecord> getSampleGenotypeRecords();
	
}
