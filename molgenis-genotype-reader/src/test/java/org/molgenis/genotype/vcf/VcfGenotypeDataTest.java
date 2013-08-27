package org.molgenis.genotype.vcf;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertNotNull;
import static org.testng.Assert.assertNull;
import static org.testng.Assert.assertTrue;

import java.io.IOException;
import java.net.URISyntaxException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import org.molgenis.genotype.Allele;

import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.ResourceTest;
import org.molgenis.genotype.Sequence;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.VcfAnnotation;
import org.molgenis.genotype.util.Utils;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variantFilter.VariantIdIncludeFilter;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

public class VcfGenotypeDataTest extends ResourceTest
{
	private VcfGenotypeData genotypeData;

	@BeforeClass
	public void beforeClass() throws IOException, URISyntaxException
	{
		genotypeData = new VcfGenotypeData(getTestVcfGz(), getTestVcfGzTbi());
	}

	@Test
	public void getSeqNames()
	{
		List<String> seqNames = genotypeData.getSeqNames();
		assertNotNull(seqNames);
		assertEquals(seqNames.size(), 3);
		assertEquals(seqNames.get(0), "1");
		assertEquals(seqNames.get(1), "2");
		assertEquals(seqNames.get(2), "3");
	}

	@Test
	public void getVariantAnnotationsMap()
	{
		Map<String, Annotation> annotations = genotypeData.getVariantAnnotationsMap();
		assertNotNull(annotations);
		assertEquals(annotations.size(), 21);

		Annotation annotation = annotations.get("NS");
		assertNotNull(annotation);
		assertEquals(annotation.getId(), "NS");
		assertTrue(annotation instanceof VcfAnnotation);
	}

	@Test
	public void testGetSequences()
	{
		List<Sequence> sequences = genotypeData.getSequences();
		assertNotNull(sequences);
		assertEquals(sequences.size(), 3);
	}

	@Test
	public void testGetSequenceByName() throws IOException
	{
		Sequence sequence = genotypeData.getSequenceByName("2");
		assertNotNull(sequence);
		assertEquals(sequence.getName(), "2");

		List<GeneticVariant> variants = Utils.iteratorToList(genotypeData.getSequenceGeneticVariants(sequence.getName()).iterator());

		assertNotNull(variants);
		assertEquals(variants.size(), 1);
		GeneticVariant variant = variants.get(0);
		assertEquals(variant.getPrimaryVariantId(), "rs4908464");
		assertEquals(variant.getStartPos(), 7569187);
		assertEquals(variant.getRefAllele().getAlleleAsString(), "G");
		assertEquals(variant.getSequenceName(), "2");
		assertEquals(variant.isSnp(), true);

		List<String> alleles = variant.getVariantAlleles().getAllelesAsString();
		assertNotNull(alleles);
		assertEquals(alleles.size(), 2);
		assertEquals(alleles.get(0), "G");
		assertEquals(alleles.get(1), "C");

		List<Alleles> sampleVariants = variant.getSampleVariants();
		assertNotNull(sampleVariants);
		assertEquals(sampleVariants.size(), 1);
		assertNotNull(sampleVariants.get(0).getAlleles());
		assertEquals(sampleVariants.get(0).getAlleles().size(), 2);
		assertEquals(sampleVariants.get(0).getAlleles().get(0).getAlleleAsString(), "C");
		assertEquals(sampleVariants.get(0).getAlleles().get(1).getAlleleAsString(), "C");
	}

	@Test
	public void testGetSequenceByNameNotExisting()
	{
		assertNull(genotypeData.getSequenceByName("Bogus"));
	}

	@Test
	public void testGetVariant()
	{
		List<GeneticVariant> variants = genotypeData.getVariantsByPos("1", 565286);
		assertNotNull(variants);
		assertEquals(variants.size(), 1);
		assertEquals(variants.get(0).getPrimaryVariantId(), "rs1578391");

		assertNotNull(genotypeData.getVariantsByPos("bogus", 8));
		assertTrue(genotypeData.getVariantsByPos("bogus", 8).isEmpty());
	}

	@Test
	public void testgetSequenceGeneticVariants() throws IOException
	{
		List<GeneticVariant> variants = Utils.iteratorToList(genotypeData.getSequenceGeneticVariants("1").iterator());
		assertEquals(variants.size(), 6);
	}

	@Test
	public void testSnpVariants()
	{
		GeneticVariant snpGeneticVariant = genotypeData.getSnpVariantByPos("1", 3172273);
		assertNotNull(snpGeneticVariant);
		assertEquals(snpGeneticVariant.isSnp(), true);
	}

	@Test
	public void testGeneticVariantAnnotations()
	{
		// NS=1;DP=4;AF=1.000;ANNOT=INT;GI=PRDM16;TI=NM_022114.3;PI=NP_071397.3
		List<GeneticVariant> variants = genotypeData.getVariantsByPos("1", 3171929);
		assertNotNull(variants);
		assertEquals(variants.size(), 1);

		GeneticVariant variant = variants.get(0);
		assertNotNull(variant.getAnnotationValues());
		assertEquals(variant.getAnnotationValues().size(), 7);

		Object annotationValue = variant.getAnnotationValues().get("NS");
		assertNotNull(annotationValue);
		assertEquals(annotationValue, Integer.valueOf(1));

		annotationValue = variant.getAnnotationValues().get("AF");
		assertNotNull(annotationValue);
		assertTrue(annotationValue instanceof List);
		@SuppressWarnings("unchecked")
		List<Float> floats = (List<Float>) annotationValue;
		assertEquals(floats.size(), 1);
		assertEquals(floats.get(0).floatValue(), 1.0, 0.001);

		annotationValue = variant.getAnnotationValues().get("ANNOT");
		assertNotNull(annotationValue);
		assertTrue(annotationValue instanceof List);
		@SuppressWarnings("unchecked")
		List<String> strings = (List<String>) annotationValue;
		assertEquals(strings.size(), 1);
		assertEquals(strings.get(0), "INT");
	}

	@Test
	public void testStopPos()
	{
		List<GeneticVariant> variants = genotypeData.getVariantsByPos("1", 565286);
		assertNotNull(variants);
		assertEquals(variants.size(), 1);

		variants = genotypeData.getVariantsByPos("3", 7569);
		assertNotNull(variants);
		assertEquals(variants.size(), 1);
	}

	@Test
	public void testGetSamplePhasing()
	{
		GeneticVariant variant = genotypeData.getVariantsByPos("1", 565286).iterator().next();
		assertEquals(genotypeData.getSamplePhasing(variant), Arrays.asList(false));
	}
	
	@Test
	public void getVariantIdMap(){
		HashMap<String, GeneticVariant> variantMap = genotypeData.getVariantIdMap();
		
		assertEquals(variantMap.size(), 7);
		
		assertTrue(variantMap.containsKey("rs35434908"));
		
		assertEquals(variantMap.get("rs35434908").getVariantAlleles().get(1), Allele.create("GTTTCA"));
		
		variantMap = genotypeData.getVariantIdMap(new VariantIdIncludeFilter("rs35434908", "rs1578391"));

		assertEquals(variantMap.size(), 2);
		
		assertTrue(variantMap.containsKey("rs35434908"));
		assertTrue(variantMap.containsKey("rs1578391"));
		
	}
	
}
