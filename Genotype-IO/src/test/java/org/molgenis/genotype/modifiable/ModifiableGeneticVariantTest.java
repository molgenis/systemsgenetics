package org.molgenis.genotype.modifiable;

import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.bgen.BgenGenotypeDataTest;
import org.molgenis.genotype.util.CalledDosageConvertor;
import org.molgenis.genotype.util.FixedSizeIterable;
import org.molgenis.genotype.util.ProbabilitiesConvertor;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.GeneticVariantMeta;
import org.molgenis.genotype.variant.GenotypeRecord;
import org.molgenis.genotype.variant.ReadOnlyGeneticVariant;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantUniqueIdProvider;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.util.*;

import static org.mockito.Mockito.mock;
import static org.testng.Assert.assertEquals;

public class ModifiableGeneticVariantTest {
    public static final ModifiableGenotypeData dummyModifiableGenotypeData = new ModifiableGenotypeDataInMemory(null);

    private GeneticVariantMeta variantMeta;

    @BeforeMethod
    public void setUp() {
        this.variantMeta = mock(GeneticVariantMeta.class);
    }

    @Test
    public void TestSampleProbabilitiesGettersWithUpdatedRef() {
		DummyVariantsProvider sampleVariantsProvider = new DummyVariantsProvider();

        Map<GeneticVariant, double[][][]> variants = new LinkedHashMap<>();
		GeneticVariant variant1 = ReadOnlyGeneticVariant.createSnp(variantMeta, "Rs1", 1, "chr1", sampleVariantsProvider, 'A', 'T', 'A');
		GeneticVariant variant2 = ReadOnlyGeneticVariant.createSnp(variantMeta, "Rs2", 20, "chr2", sampleVariantsProvider, 'G', 'T', 'G');

		variants.put(variant1, new double[][][]{{{0, 1}, {0, 1}}, {{1, 0}, {1, 0}}, {{0, 1}, {1, 0}}});
		variants.put(variant2, new double[][][]{
				{{0, 1}, {0.6, 0.4}, {0.5, 0.5}, {1, 0}},
				{{0.5, 0.5}, {0.5, 0.5}, {0.2, 0.8}, {0, 1}},
				{{0.5, 0.5}, {0.5, 0.5}, {0.8, 0.2}, {1, 0}}});

		sampleVariantsProvider.setVariants(variants);

        HashSet<ModifiableGeneticVariant> excludeList = new HashSet<>();

		List<ModifiableGeneticVariant> modifiableVariants = new ArrayList<>();
        Iterable<ModifiableGeneticVariant> modifiableVariantsIterator = ModifiableGeneticVariantIterator
                .createModifiableGeneticVariantIterable(variants.keySet().iterator(), dummyModifiableGenotypeData, excludeList);

		modifiableVariantsIterator.forEach(modifiableVariants::add);

		assertEquals(modifiableVariants.get(0).getOriginalVariant(), variant1);
		assertEquals(modifiableVariants.get(1).getOriginalVariant(), variant2);

        dummyModifiableGenotypeData.updateRefAllele(modifiableVariants.get(0), Allele.T);
        dummyModifiableGenotypeData.updateRefAllele(modifiableVariants.get(1), Allele.T);

		double[][] probsByModifier = modifiableVariants.get(0).getSampleGenotypeProbabilitiesComplex();
		double[][] probsByModifier1 = modifiableVariants.get(1).getSampleGenotypeProbabilitiesComplex();

		assertEquals(probsByModifier, new double[][]{{1,0,0}, {0,0,1}, {0,1,0}});
		BgenGenotypeDataTest.assertProbabilityEquality(probsByModifier1, new double[][]{
				{0, 0.2, 0.5, 0.3, 0},
				{0.2, 0.45, 0.3, 0.05, 0.0},
				{0.0, 0.05, 0.3, 0.45, 0.2}}, 1E-7);

		double[][][] sampleGenotypeProbabilitiesPhased = modifiableVariants.get(0).getSampleGenotypeProbabilitiesPhased();

		assertEquals(sampleGenotypeProbabilitiesPhased,
				new double[][][]{{{1, 0}, {1, 0}}, {{0, 1}, {0, 1}}, {{1, 0}, {0, 1}}});
		assertEquals(
				modifiableVariants.get(1).getSampleGenotypeProbabilitiesPhased(),
				new double[][][]{
						{{1, 0}, {0.4, 0.6}, {0.5, 0.5}, {0, 1}},
						{{0.5, 0.5}, {0.5, 0.5}, {0.8, 0.2}, {1, 0}},
						{{0.5, 0.5}, {0.5, 0.5}, {0.2, 0.8}, {0, 1}}}
				);
    }

    private static class DummyVariantsProvider implements SampleVariantsProvider {

        private Map<GeneticVariant, double[][][]> variants;
        private final int sampleVariantProviderUniqueId;
		private static final double DEFAULT_MINIMUM_POSTERIOR_PROBABILITY_TO_CALL = 0.4f;

		public DummyVariantsProvider() {
            super();
			sampleVariantProviderUniqueId = SampleVariantUniqueIdProvider.getNextUniqueId();
        }

        @Override
        public List<Alleles> getSampleVariants(GeneticVariant variant) {
            return ProbabilitiesConvertor.convertProbabilitiesToAlleles(
					variant.getSampleGenotypeProbilities(),
					variant.getVariantAlleles(),
					DEFAULT_MINIMUM_POSTERIOR_PROBABILITY_TO_CALL);
        }

        @Override
        public int cacheSize() {
            return 0;
        }

        @Override
        public List<Boolean> getSamplePhasing(GeneticVariant variant) {
			Boolean[] array = new Boolean[variant.getSampleGenotypeProbabilitiesComplex().length];
			Arrays.fill(array, Boolean.TRUE);
            return Arrays.asList(array);
        }

        @Override
        public int getSampleVariantProviderUniqueId() {
            return sampleVariantProviderUniqueId;
        }

        @Override
        public byte[] getSampleCalledDosage(GeneticVariant variant) {
            return CalledDosageConvertor.convertCalledAllelesToCalledDosage(getSampleVariants(variant),
                    variant.getVariantAlleles(), variant.getRefAllele());
        }

        @Override
        public float[] getSampleDosage(GeneticVariant variant) {
            return CalledDosageConvertor.convertCalledAllelesToDosage(getSampleVariants(variant),
                    variant.getVariantAlleles(), variant.getRefAllele());
        }

        @Override
        public float[][] getSampleProbilities(GeneticVariant variant) {
            return ProbabilitiesConvertor.convertBiallelicComplexProbabilitiesToProbabilities(
            		getSampleProbabilitiesComplex(variant));
        }

        @Override
        public double[][] getSampleProbabilitiesComplex(GeneticVariant variant) {
            return ProbabilitiesConvertor.convertPhasedProbabilitiesToComplexProbabilities(
					getSampleProbabilitiesPhased(variant));
        }

        @Override
        public double[][][] getSampleProbabilitiesPhased(GeneticVariant variant) {
            return variants.get(variant);
        }

        @Override
        public FixedSizeIterable<GenotypeRecord> getSampleGenotypeRecords(GeneticVariant variant) {
            return null;
        }

		public void setVariants(Map<GeneticVariant, double[][][]> variants) {
			this.variants = variants;
		}
	}
}

