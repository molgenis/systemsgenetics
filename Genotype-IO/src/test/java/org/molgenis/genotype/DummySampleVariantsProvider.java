package org.molgenis.genotype;

import java.util.List;

import org.molgenis.genotype.util.CalledDosageConvertor;
import org.molgenis.genotype.util.FixedSizeIterable;
import org.molgenis.genotype.util.ProbabilitiesConvertor;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.GenotypeRecord;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantUniqueIdProvider;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;

public class DummySampleVariantsProvider implements SampleVariantsProvider
{

	private final List<Alleles> variantAlleles;
	private final int sampleVariantProviderUniqueId;

	public DummySampleVariantsProvider(List<Alleles> variantAlleles)
	{
		super();
		this.variantAlleles = variantAlleles;
		sampleVariantProviderUniqueId = SampleVariantUniqueIdProvider.getNextUniqueId();
	}

	@Override
	public List<Alleles> getSampleVariants(GeneticVariant variant)
	{
		return variantAlleles;
	}

	@Override
	public int cacheSize()
	{
		return 0;
	}

	@Override
	public List<Boolean> getSamplePhasing(GeneticVariant variant)
	{
		return null;
	}

	@Override
	public boolean arePhasedProbabilitiesPresent(GeneticVariant variant) {
		return false;
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
	
	@Override
	public float[][] getSampleProbilities(GeneticVariant variant) {
		return ProbabilitiesConvertor.convertCalledAllelesToProbability(variant.getSampleVariants(), variant.getVariantAlleles());
	}

	@Override
	public double[][] getSampleProbabilitiesComplex(GeneticVariant variant) {
		return ProbabilitiesConvertor.convertProbabilitiesToComplexProbabilities(getSampleProbilities(variant));
	}

	@Override
	public double[][][] getSampleProbabilitiesPhased(GeneticVariant variant) {
		throw new GenotypeDataException("Phased data not available");
	}

	@Override
	public FixedSizeIterable<GenotypeRecord> getSampleGenotypeRecords(GeneticVariant variant)
	{
		return null;
	}
}
