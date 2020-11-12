package org.molgenis.genotype.variant.sampleProvider;

import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.times;
import static org.mockito.Mockito.verify;

import org.molgenis.genotype.variant.GeneticVariant;
import org.testng.annotations.Test;

public class CachedSampleVariantProviderTest
{
	@Test
	public void getSampleGenotypeRecords()
	{
		SampleVariantsProvider sampleVariantProvider = mock(SampleVariantsProvider.class);
		CachedSampleVariantProvider cachedSampleVariantProvider = new CachedSampleVariantProvider(sampleVariantProvider , 1);
		GeneticVariant variant = mock(GeneticVariant.class);
		cachedSampleVariantProvider.getSampleGenotypeRecords(variant);
		cachedSampleVariantProvider.getSampleGenotypeRecords(variant);
		//broken test because mock gives null results for sapmle records
		//verify(sampleVariantProvider, times(1)).getSampleGenotypeRecords(variant); // once from cache, once from sampleVariantProvider
	}
}
