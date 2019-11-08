package org.molgenis.genotype.variant.sampleProvider;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.util.FixedSizeIterable;
import org.molgenis.genotype.util.ProbabilitiesConvertor;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.GenotypeRecord;

public class SwappingSampleVariantsProvider implements SampleVariantsProvider
{
	private final SampleVariantsProvider sampleVariantsProvider;
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

	@Override
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

	@Override
	public float[][] getSampleProbilities(GeneticVariant variant) {
		return sampleVariantsProvider.getSampleProbilities(variant);
	}

	@Override
	public double[][] getSampleGenotypeProbabilitiesBgen(GeneticVariant variant) {
		return ProbabilitiesConvertor.convertProbabilitiesToBgenProbabilities(getSampleProbilities(variant));
	}

	@Override
	public double[][][] getSampleGenotypeProbabilitiesBgenPhased(GeneticVariant variant) {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	@Override
	public FixedSizeIterable<GenotypeRecord> getSampleGenotypeRecords(GeneticVariant variant)
	{
		return sampleVariantsProvider.getSampleGenotypeRecords(variant);
	}

}
