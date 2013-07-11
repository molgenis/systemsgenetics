package org.molgenis.genotype.impute2;

import java.util.ArrayList;
import java.util.List;

import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.util.CalledDosageConvertor;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantUniqueIdProvider;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;

/**
 * SampleVariantsProvider that returns the sample alleles of a HapsEntry
 * 
 * @author erwin
 * 
 */
public class HapsEntrySampleVariantsProvider implements SampleVariantsProvider
{
	private HapsEntry hapsEntry;
	private final int sampleVariantProviderUniqueId;

	public HapsEntrySampleVariantsProvider(HapsEntry hapsEntry)
	{
		this.hapsEntry = hapsEntry;
		sampleVariantProviderUniqueId = SampleVariantUniqueIdProvider.getNextUniqueId();
	}

	@Override
	public List<Alleles> getSampleVariants(GeneticVariant variant)
	{
		if (!variant.getPrimaryVariantId().equals(hapsEntry.getSnpId()))
		{
			throw new IllegalArgumentException("Ids don't match");
		}

		List<Alleles> allelesList = new ArrayList<Alleles>();
		for (String[] sampleAlleles : hapsEntry.getSampleAlleles())
		{
			Allele allele1 = createAllele(sampleAlleles[0]);
			Allele allele2 = createAllele(sampleAlleles[1]);
			Alleles alleles = Alleles.createAlleles(allele1, allele2);
			allelesList.add(alleles);
		}

		return allelesList;
	}

	/*
	 * Sample can have an asterisk directly after the allele indicating that it
	 * is unphased
	 */
	private Allele createAllele(String sample)
	{
		char allele = sample.charAt(0);
		switch (allele)
		{
			case '?':
				return Allele.ZERO;
			case '0':
				return Allele.create(hapsEntry.getFirstAllele());
			case '1':
				return Allele.create(hapsEntry.getSecondAllele());
			default:
				throw new IllegalArgumentException("[" + sample + "] is an invalid value for a haps sample value");
		}

	}

	/*
	 * Sample can have an asterisk directly after the allele indicating that it
	 * is unphased
	 */
	@Override
	public List<Boolean> getSamplePhasing(GeneticVariant variant)
	{
		if (!variant.getPrimaryVariantId().equals(hapsEntry.getSnpId()))
		{
			throw new IllegalArgumentException("Ids don't match");
		}

		List<Boolean> phasingList = new ArrayList<Boolean>();
		for (String[] sampleAlleles : hapsEntry.getSampleAlleles())
		{
			boolean phased = !sampleAlleles[0].endsWith("*");
			phasingList.add(phased);
		}

		return phasingList;
	}

	@Override
	public int cacheSize()
	{
		return 0;
	}

	@Override
	public int getSampleVariantProviderUniqueId()
	{
		return sampleVariantProviderUniqueId;
	}

	@Override
	public byte[] getSampleCalledDosage(GeneticVariant variant)
	{
		return CalledDosageConvertor.convertCalledAllelesToCalledDosage(getSampleVariants(variant),
				variant.getVariantAlleles(), variant.getRefAllele());
	}

	@Override
	public float[] getSampleDosage(GeneticVariant variant)
	{
		return CalledDosageConvertor.convertCalledAllelesToDosage(getSampleVariants(variant),
				variant.getVariantAlleles(), variant.getRefAllele());
	}

}
