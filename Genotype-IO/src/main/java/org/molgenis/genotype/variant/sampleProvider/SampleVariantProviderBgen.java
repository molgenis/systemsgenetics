package org.molgenis.genotype.variant.sampleProvider;

import org.molgenis.genotype.variant.ReadOnlyGeneticVariantBgen;

public interface SampleVariantProviderBgen extends SampleVariantsProvider {
    ReadOnlyGeneticVariantBgen extendReadOnlyGeneticVariantBgen(ReadOnlyGeneticVariantBgen variant);
}
