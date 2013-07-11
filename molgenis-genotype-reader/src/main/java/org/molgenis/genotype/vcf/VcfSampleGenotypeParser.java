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
public class VcfSampleGenotypeParser
{
	private static final char PHASED_SEPARATOR = '|';
	private static final char UNPHASED_SEPARATOR = '/';
	private String genotype;

	public VcfSampleGenotypeParser(String genotype)
	{
		if (genotype == null) throw new IllegalArgumentException("Genotype is null");
		this.genotype = genotype;
	}

	public VcfSampleGenotype parse()
	{
		//Length will almost always be 2.
		List<Character> alleleIndices = new ArrayList<Character>(2);
		
		List<Boolean> phasing = new ArrayList<Boolean>();

		for (char c : genotype.toCharArray())
		{
			switch (c)
			{
				case PHASED_SEPARATOR:
					phasing.add(true);
					break;
				case UNPHASED_SEPARATOR:
					phasing.add(false);
					break;
				default:
					alleleIndices.add(c);
			}

		}

		return new VcfSampleGenotype(alleleIndices, phasing);
	}
}
