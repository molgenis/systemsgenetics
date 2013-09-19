package org.molgenis.genotype.variant.sampleProvider;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.variant.GeneticVariant;

public class SwappingSampleVariantsProvider implements SampleVariantsProvider
{
	private SampleVariantsProvider sampleVariantsProvider;
	private final int sampleVariantProviderUniqueId;

	public SwappingSampleVariantsProvider(SampleVariantsProvider sampleVariantsProvider)
	{
		this.sampleVariantsProvider = sampleVariantsProvider;
		sampleVariantProviderUniqueId = SampleVariantUniqueIdProvider.getNextUniqueId();
	}

	@Override
	public List<Alleles> getSampleVariants(GeneticVariant variant)
	{

		List<Alleles> sampleAlleles = sampleVariantsProvider.getSampleVariants(variant);
		List<Alleles> swapped = new ArrayList<Alleles>(sampleAlleles.size());
		for (Alleles alleles : sampleAlleles)
		{
			swapped.add(alleles.getComplement());
		}

		return Collections.unmodifiableList(swapped);
	}

	@Override
	public int cacheSize()
	{
		return 0;
	}

	@Override
	public List<Boolean> getSamplePhasing(GeneticVariant variant)
	{
		return sampleVariantsProvider.getSamplePhasing(variant);
	}

	public int getSampleVariantProviderUniqueId()
	{
		return sampleVariantProviderUniqueId;
	}

	@Override
	public byte[] getSampleCalledDosage(GeneticVariant variant)
	{
		return sampleVariantsProvider.getSampleCalledDosage(variant);
	}

	@Override
	public float[] getSampleDosage(GeneticVariant variant)
	{
		return sampleVariantsProvider.getSampleDosage(variant);
	}

}
