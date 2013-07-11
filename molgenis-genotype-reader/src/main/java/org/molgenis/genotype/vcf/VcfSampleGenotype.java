package org.molgenis.genotype.vcf;

import java.util.ArrayList;
import java.util.List;

/**
 * GT : genotype, encoded as allele values separated by either of ”/” or “|”.
 * The allele values are 0 for the reference allele (what is in the REF field),
 * 1 for the first allele listed in ALT, 2 for the second allele list in ALT and
 * so on. For diploid calls examples could be 0/1, 1|0, or 1/2, etc. For haploid
 * calls, e.g. on Y, male non-pseudoautosomal X, or mitochondrion, only one
 * allele value should be given; a triploid call might look like 0/0/1. If a
 * call cannot be made for a sample at a given locus, ”.”should be specified for
 * each missing allele in the GT field (for example "./." for a diploid genotype
 * and "." for haploid genotype). The meanings of the separators are as follows
 * (see the PS field below for more details on incorporating phasing information
 * into the genotypes):
 * 
 * / : genotype unphased | : genotype phased
 * 
 * 
 * 
 * @author erwin
 * 
 */
public class VcfSampleGenotype
{
	private List<Character> alleleIndices;
	private List<Boolean> phasing;

	public VcfSampleGenotype(List<Character> alleleIndices, List<Boolean> phasing)
	{
		if (alleleIndices == null) throw new IllegalArgumentException("alleleIndices list is null");
		if (phasing == null) throw new IllegalArgumentException("Phasing list is null");

		this.alleleIndices = alleleIndices;
		this.phasing = phasing;
	}

	/**
	 * Get the allele values '.' is used for unknown values
	 * 
	 * @return
	 */
	public List<Character> getAlleleIndices()
	{
		return alleleIndices;
	}

	public int getAlleleIndexCount()
	{
		return alleleIndices.size();
	}

	public List<Boolean> getPhasing()
	{
		return phasing;
	}

	/**
	 * Get all sample variants (alleles). Can contain null values !!!! if
	 * unknown
	 */
	public List<String> getSamleVariants(List<String> alleles)
	{
		List<String> sampleVariants = new ArrayList<String>(alleleIndices.size());

		for (Character index : alleleIndices)
		{
			String variant = null;
			String value = VcfUtils.checkNullValue(Character.toString(index));
			if (value != null)
			{
				variant = alleles.get(Integer.parseInt(value));
			}

			sampleVariants.add(variant);
		}

		return sampleVariants;
	}
}
