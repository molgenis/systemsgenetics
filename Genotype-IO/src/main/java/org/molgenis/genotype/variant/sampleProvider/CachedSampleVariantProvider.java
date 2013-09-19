package org.molgenis.genotype.variant.sampleProvider;

import java.util.List;

import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.util.Cache;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 * Cached sample variant provider to prevent reloading a SNPs that is accessed
 * multiple times in a sort periode.
 * 
 * @author Patrick Deelen
 * 
 */
public class CachedSampleVariantProvider implements SampleVariantsProvider
{

	private final SampleVariantsProvider sampleVariantProvider;
	private final Cache<GeneticVariant, List<Alleles>> cache;
	private final Cache<GeneticVariant, List<Boolean>> phasingCache;
	private final Cache<GeneticVariant, byte[]> calledDosageCache;
	private final Cache<GeneticVariant, float[]> dosageCache;
	private final int cacheSize;
	private final int sampleVariantProviderUniqueId;

	public CachedSampleVariantProvider(SampleVariantsProvider sampleVariantProvider, int cacheSize)
	{
		this.sampleVariantProvider = sampleVariantProvider;
		this.cache = new Cache<GeneticVariant, List<Alleles>>(cacheSize);
		this.phasingCache = new Cache<GeneticVariant, List<Boolean>>(cacheSize);
		this.calledDosageCache = new Cache<GeneticVariant, byte[]>(cacheSize);
		this.dosageCache = new Cache<GeneticVariant, float[]>(cacheSize);
		this.cacheSize = cacheSize;
		sampleVariantProviderUniqueId = SampleVariantUniqueIdProvider.getNextUniqueId();
	}

	@Override
	public List<Alleles> getSampleVariants(GeneticVariant variant)
	{
		if (cache.containsKey(variant))
		{
			return cache.get(variant);
		}
		else
		{
			List<Alleles> variantAlleles = sampleVariantProvider.getSampleVariants(variant);
			cache.put(variant, variantAlleles);
			return variantAlleles;
		}

	}

	@Override
	public int cacheSize()
	{
		return cacheSize;
	}

	@Override
	public List<Boolean> getSamplePhasing(GeneticVariant variant)
	{
		if (phasingCache.containsKey(variant))
		{
			return phasingCache.get(variant);
		}

		List<Boolean> phasing = sampleVariantProvider.getSamplePhasing(variant);
		phasingCache.put(variant, phasing);
		return phasing;
	}

	@Override
	public int getSampleVariantProviderUniqueId()
	{
		return sampleVariantProviderUniqueId;
	}

	@Override
	public byte[] getSampleCalledDosage(GeneticVariant variant)
	{
		if (calledDosageCache.containsKey(variant))
		{
			return calledDosageCache.get(variant);
		}

		byte[] calledDosage = sampleVariantProvider.getSampleCalledDosage(variant);
		calledDosageCache.put(variant, calledDosage);
		return calledDosage;
	}

	@Override
	public float[] getSampleDosage(GeneticVariant variant)
	{
		if (dosageCache.containsKey(variant))
		{
			return dosageCache.get(variant);
		}

		float[] dosage = sampleVariantProvider.getSampleDosage(variant);
		dosageCache.put(variant, dosage);
		return dosage;
	}
}
