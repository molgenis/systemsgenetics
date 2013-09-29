package org.molgenis.genotype.vcf;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertNotNull;
import static org.testng.Assert.assertNull;
import static org.testng.Assert.assertTrue;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.probabilities.SampleVariantProbabilities;
import org.molgenis.genotype.util.CalledDosageConvertor;
import org.molgenis.genotype.util.ProbabilitiesConvertor;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.VariantLineMapper;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantUniqueIdProvider;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;
import org.testng.annotations.Test;

public class VcfVariantLineMapperTest implements SampleVariantsProvider {

	private static final List<String> COL_NAMES = Arrays.asList("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
			"INFO", "FORMAT", "sample1");
	private final int sampleVariantProviderUniqueId;

	public VcfVariantLineMapperTest() {
		sampleVariantProviderUniqueId = SampleVariantUniqueIdProvider.getNextUniqueId();
	}

	@Test
	public void mapLineSnp() {
		String line = "1	565286	rs1578391	C	T	.	flt	NS=1;DP=5;AF=1.000;ANNOT=INT;GI=LOC100131754	GT:DP:EC:CONFS	1/1:5:5:5.300,5.300,1.000,1.000,1.000,1.000,1.000";

		VariantLineMapper mapper = new VcfVariantLineMapper(COL_NAMES, Collections.<Annotation>emptyList(),
				Collections.<String, String>emptyMap(), this);
		GeneticVariant variant = mapper.mapLine(line);
		assertNotNull(variant);
		assertEquals(variant.isSnp(), true);
		assertEquals(variant.getSequenceName(), "1");
		assertEquals(variant.getStartPos(), 565286);
		assertEquals(variant.getPrimaryVariantId(), "rs1578391");

		List<String> ids = variant.getAllIds();
		assertNotNull(ids);
		assertEquals(ids.size(), 1);
		assertEquals(ids.get(0), "rs1578391");

		assertEquals(variant.getRefAllele().getAlleleAsString(), "C");

		List<String> alleles = variant.getVariantAlleles().getAllelesAsString();
		assertNotNull(alleles);
		assertEquals(alleles.size(), 2);
		assertEquals(alleles.get(0), "C");
		assertEquals(alleles.get(1), "T");

		assertEquals(variant.getRefAllele().getAlleleAsString(), "C");
		List<String> snpAlleles = variant.getVariantAlleles().getAllelesAsString();
		assertNotNull(snpAlleles);
		assertEquals(snpAlleles.size(), 2);
		assertEquals(snpAlleles.get(0), "C");
		assertEquals(snpAlleles.get(1), "T");
	}

	@Test
	public void mapLineInsert() {
		String line = "1	565286	.	C	CTA,CA	.	flt	NS=1;DP=5;AF=1.000;ANNOT=INT;GI=LOC100131754	GT:DP:EC:CONFS	1/1:5:5:5.300,5.300,1.000,1.000,1.000,1.000,1.000";

		VariantLineMapper mapper = new VcfVariantLineMapper(COL_NAMES, Collections.<Annotation>emptyList(),
				Collections.<String, String>emptyMap(), this);
		GeneticVariant variant = mapper.mapLine(line);
		assertNotNull(variant);

		assertEquals(variant.getSequenceName(), "1");
		assertEquals(variant.getStartPos(), 565286);
		assertNull(variant.getPrimaryVariantId());

		List<String> ids = variant.getAlternativeVariantIds();
		assertNotNull(ids);
		assertTrue(ids.isEmpty());

		assertEquals(variant.getRefAllele().getAlleleAsString(), "C");

		List<String> alleles = variant.getVariantAlleles().getAllelesAsString();
		assertNotNull(alleles);
		assertEquals(alleles.size(), 3);
		assertTrue(alleles.contains("C"));
		assertTrue(alleles.contains("CTA"));
		assertTrue(alleles.contains("CA"));
	}

	@Override
	public List<Alleles> getSampleVariants(GeneticVariant variant) {
		return Collections.emptyList();
	}

	@Override
	public int cacheSize() {
		return 0;
	}

	@Override
	public List<Boolean> getSamplePhasing(GeneticVariant variant) {
		return null;
	}

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
	public SampleVariantProbabilities[] getSampleProbilities(GeneticVariant variant) {
		return ProbabilitiesConvertor.convertCalledAllelesToProbability(variant.getSampleVariants(), variant.getVariantAlleles());
	}
}
