package org.molgenis.genotype.variant.sampleProvider;

import java.util.List;

import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.util.Cache;
import org.molgenis.genotype.util.FixedSizeIterable;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.GenotypeRecord;

/**
 * Cached sample variant provider to prevent reloading a SNPs that is accessed
 * multiple times in a short period.
 *
 * @author Patrick Deelen
 *
 */
public class CachedSampleVariantProvider implements SampleVariantsProvider {

	private final SampleVariantsProvider sampleVariantProvider;
	private final Cache<GeneticVariant, List<Alleles>> cache;
	private final Cache<GeneticVariant, List<Boolean>> phasingCache;
	private final Cache<GeneticVariant, Boolean> phasedProbCache;
	private final Cache<GeneticVariant, byte[]> calledDosageCache;
	private final Cache<GeneticVariant, float[]> dosageCache;
	private final Cache<GeneticVariant, float[][]> probCache;
	private final Cache<GeneticVariant, double[][]> probCacheComplex;
	private final Cache<GeneticVariant, double[][][]> probCachePhased;
	private final Cache<GeneticVariant, FixedSizeIterable<GenotypeRecord>> genotypeRecordCache;
	private final int cacheSize;
	private final int sampleVariantProviderUniqueId;

	public CachedSampleVariantProvider(SampleVariantsProvider sampleVariantProvider, int cacheSize) {
		this.sampleVariantProvider = sampleVariantProvider;
		this.cache = new Cache<GeneticVariant, List<Alleles>>(cacheSize);
		this.phasingCache = new Cache<GeneticVariant, List<Boolean>>(cacheSize);
		this.phasedProbCache = new Cache<GeneticVariant, Boolean>(cacheSize);
		this.calledDosageCache = new Cache<GeneticVariant, byte[]>(cacheSize);
		this.dosageCache = new Cache<GeneticVariant, float[]>(cacheSize);
		this.probCache = new Cache<GeneticVariant, float[][]>(cacheSize);
		this.probCacheComplex = new Cache<GeneticVariant, double[][]>(cacheSize);
		this.probCachePhased = new Cache<GeneticVariant, double[][][]>(cacheSize);
		this.genotypeRecordCache = new Cache<GeneticVariant, FixedSizeIterable<GenotypeRecord>>(cacheSize);
		this.cacheSize = cacheSize;
		sampleVariantProviderUniqueId = SampleVariantUniqueIdProvider.getNextUniqueId();
	}

	@Override
	public List<Alleles> getSampleVariants(GeneticVariant variant) {
		List<Alleles> variantAlleles = cache.get(variant);
		if (variantAlleles == null) {
			variantAlleles = sampleVariantProvider.getSampleVariants(variant);
			cache.put(variant, variantAlleles);
		}
		return variantAlleles;

	}

	@Override
	public int cacheSize() {
		return cacheSize;
	}

	@Override
	public List<Boolean> getSamplePhasing(GeneticVariant variant) {

		List<Boolean> phasing = phasingCache.get(variant);
		if (phasing == null) {
			phasing = sampleVariantProvider.getSamplePhasing(variant);
			phasingCache.put(variant, phasing);
		}
		return phasing;
	}

	@Override

	public boolean arePhasedProbabilitiesPresent(GeneticVariant variant) {
		if (phasedProbCache.containsKey(variant))
		{
			return phasedProbCache.get(variant);
		}

		boolean phasing = sampleVariantProvider.arePhasedProbabilitiesPresent(variant);
		phasedProbCache.put(variant, phasing);
		return phasing;
	}

	@Override
	public int getSampleVariantProviderUniqueId()
	{
		return sampleVariantProviderUniqueId;
	}

	@Override
	public byte[] getSampleCalledDosage(GeneticVariant variant) {

		byte[] calledDosage = calledDosageCache.get(variant);
		if (calledDosage == null) {
			calledDosage = sampleVariantProvider.getSampleCalledDosage(variant);
			calledDosageCache.put(variant, calledDosage);
		}
		return calledDosage;
	}

	@Override
	public float[] getSampleDosage(GeneticVariant variant) {

		float[] dosage = dosageCache.get(variant);
		if (dosage == null) {
			dosage = sampleVariantProvider.getSampleDosage(variant);
			dosageCache.put(variant, dosage);
		}

		return dosage;
	}

	@Override
	public float[][] getSampleProbilities(GeneticVariant variant) {

		float[][] probs = probCache.get(variant);
		if (probs == null) {
			probs = sampleVariantProvider.getSampleProbilities(variant);
			probCache.put(variant, probs);
		}
		return probs;
	}

	@Override
	public double[][] getSampleProbabilitiesComplex(GeneticVariant variant) {

		double[][] probs = probCacheComplex.get(variant);
		if (probs == null) {
			probs = sampleVariantProvider.getSampleProbabilitiesComplex(variant);
			probCacheComplex.put(variant, probs);
		}
		return probs;

	}

	@Override
	public double[][][] getSampleProbabilitiesPhased(GeneticVariant variant) {

		double[][][] probs = probCachePhased.get(variant);
		if (probs == null) {
			probs = sampleVariantProvider.getSampleProbabilitiesPhased(variant);
			probCachePhased.put(variant, probs);
		}
		return probs;
	}

	@Override
	public FixedSizeIterable<GenotypeRecord> getSampleGenotypeRecords(GeneticVariant variant) {

		FixedSizeIterable<GenotypeRecord> sampleGenotypeRecords = genotypeRecordCache.get(variant);
		if (sampleGenotypeRecords == null) {
			sampleGenotypeRecords = sampleVariantProvider.getSampleGenotypeRecords(variant);
			genotypeRecordCache.put(variant, sampleGenotypeRecords);
		}

		return sampleGenotypeRecords;
	}
}
