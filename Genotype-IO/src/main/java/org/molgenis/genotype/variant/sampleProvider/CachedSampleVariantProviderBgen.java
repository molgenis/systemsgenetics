package org.molgenis.genotype.variant.sampleProvider;

import org.molgenis.genotype.variant.ReadOnlyGeneticVariantBgen;

public class CachedSampleVariantProviderBgen extends CachedSampleVariantProvider implements SampleVariantProviderBgen {
    protected final SampleVariantProviderBgen sampleVariantProvider;

    public CachedSampleVariantProviderBgen(SampleVariantProviderBgen sampleVariantProvider, int cacheSize) {
        super(sampleVariantProvider, cacheSize);
        this.sampleVariantProvider = sampleVariantProvider;
    }

    @Override
    public ReadOnlyGeneticVariantBgen extendReadOnlyGeneticVariantBgen(ReadOnlyGeneticVariantBgen variant) {
        return sampleVariantProvider.extendReadOnlyGeneticVariantBgen(variant);
    }
}
