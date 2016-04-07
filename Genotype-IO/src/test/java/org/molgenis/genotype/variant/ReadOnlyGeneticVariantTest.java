package org.molgenis.genotype.variant;

import static org.mockito.Mockito.mock;
import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertNotNull;
import static org.testng.Assert.assertTrue;
import static org.testng.AssertJUnit.assertNull;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.DummySampleVariantsProvider;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

public class ReadOnlyGeneticVariantTest
{
	private GeneticVariant testInstance;
	private GeneticVariantMeta variantMeta;
	
	@BeforeClass
	public void setUp()
	{

		ArrayList<Alleles> sampleAlleles = new ArrayList<Alleles>();
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('C', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		SampleVariantsProvider sampleAllelesProvider = new DummySampleVariantsProvider(sampleAlleles);

		variantMeta = mock(GeneticVariantMeta.class);
		testInstance = ReadOnlyGeneticVariant.createSnp(variantMeta, "rs1", 1, "chr1", sampleAllelesProvider, 'A', 'C');
	}

	@Test
	public void getAllIds()
	{
		ArrayList<String> expected = new ArrayList<String>();
		expected.add("rs1");
		assertEquals(testInstance.getAllIds(), expected);
	}

	@Test
	public void getAlleleCount()
	{
		assertEquals(testInstance.getAlleleCount(), 2);
	}

	@Test
	public void getCalledDosages()
	{
		byte[] expected =
		{ 1, 1, 0, 2, 2 };
		assertEquals(testInstance.getSampleCalledDosages(), expected);
	}

	@Test
	public void getMinorAllele()
	{
		assertEquals(testInstance.getMinorAllele(), Allele.C);
	}

	@Test
	public void getMinorAlleleFrequency()
	{
		assertEquals(testInstance.getMinorAlleleFrequency(), 0.4d, 0.00001);
	}

	@Test
	public void getRefAllele()
	{
		assertNull(testInstance.getRefAllele());
	}
	
	@Test
	public void getAlternativeAlleles()
	{
		assertEquals(testInstance.getAlternativeAlleles(), Alleles.createBasedOnChars('A', 'C'));
	}

	@Test
	public void getSequenceName()
	{
		assertEquals(testInstance.getSequenceName(), "chr1");
	}

	@Test
	public void getStartPos()
	{
		assertEquals(testInstance.getStartPos(), 1);
	}

	@Test
	public void isBiallelic()
	{
		assertEquals(testInstance.isBiallelic(), true);
	}

	@Test
	public void isSnp()
	{
		assertEquals(testInstance.isSnp(), true);
	}

	@Test
	public void getPrimaryVatiantId()
	{
		GeneticVariant variant = createGeneticVariant(Arrays.asList("X"));
		assertNotNull(variant.getPrimaryVariantId());
		assertEquals("X", variant.getPrimaryVariantId());

		variant = createGeneticVariant(Arrays.asList("X", "Y"));
		assertNotNull(variant.getPrimaryVariantId());
		assertEquals("X", variant.getPrimaryVariantId());

		variant = createGeneticVariant(Arrays.asList("X", "Y", "Z"));
		assertNotNull(variant.getPrimaryVariantId());
		assertEquals("X", variant.getPrimaryVariantId());
	}

	@Test
	public void getAlternativeVariantIds()
	{
		GeneticVariant variant = createGeneticVariant(Arrays.asList("X"));
		assertNotNull(variant.getAlternativeVariantIds());
		assertTrue(variant.getAlternativeVariantIds().isEmpty());

		variant = createGeneticVariant(Arrays.asList("X", "Y"));
		assertNotNull(variant.getAlternativeVariantIds());
		assertEquals(variant.getAlternativeVariantIds().size(), 1);
		assertEquals(variant.getAlternativeVariantIds().get(0), "Y");

		variant = createGeneticVariant(Arrays.asList("X", "Y", "Z"));
		assertNotNull(variant.getAlternativeVariantIds());
		assertEquals(variant.getAlternativeVariantIds().size(), 2);
		assertEquals(variant.getAlternativeVariantIds().get(0), "Y");
		assertEquals(variant.getAlternativeVariantIds().get(1), "Z");
	}

	private GeneticVariant createGeneticVariant(List<String> ids)
	{
		return ReadOnlyGeneticVariant.createSnp(variantMeta, ids, 1, "chr1", null, 'A', 'C');
	}

	@Test
	public void isAtOrGcSnp()
	{
		assertEquals(testInstance.isAtOrGcSnp(), false);

		GeneticVariant testInstance2 = ReadOnlyGeneticVariant.createSnp(variantMeta, "rs1", 1, "chr1", null, 'G', 'C');
		GeneticVariant testInstance3 = ReadOnlyGeneticVariant.createSnp(variantMeta, "rs1", 1, "chr1", null, 'C', 'G');
		GeneticVariant testInstance4 = ReadOnlyGeneticVariant.createSnp(variantMeta, "rs1", 1, "chr1", null, 'A', 'T');
		GeneticVariant testInstance5 = ReadOnlyGeneticVariant.createSnp(variantMeta, "rs1", 1, "chr1", null, 'T', 'A');
		GeneticVariant testInstance6 = ReadOnlyGeneticVariant.createSnp(variantMeta, "rs1", 1, "chr1", null, 'C', 'A');
		GeneticVariant testInstance7 = ReadOnlyGeneticVariant.createSnp(variantMeta, "rs1", 1, "chr1", null, 'G', 'A');
		GeneticVariant testInstance8 = ReadOnlyGeneticVariant.createSnp(variantMeta, "rs1", 1, "chr1", null, 'G', 'T');

		assertEquals(testInstance2.isAtOrGcSnp(), true);
		assertEquals(testInstance3.isAtOrGcSnp(), true);
		assertEquals(testInstance4.isAtOrGcSnp(), true);
		assertEquals(testInstance5.isAtOrGcSnp(), true);
		assertEquals(testInstance6.isAtOrGcSnp(), false);
		assertEquals(testInstance7.isAtOrGcSnp(), false);
		assertEquals(testInstance8.isAtOrGcSnp(), false);

	}

	@Test(expectedExceptions = GenotypeDataException.class)
	public void checkRefIsPartOfAlleles()
	{
		@SuppressWarnings("unused")
		GeneticVariant testInstance2 = ReadOnlyGeneticVariant.createSnp(variantMeta, "rs1", 1, "chr1", null, 'A', 'C', 'T');
	}

	public void testMoveRefToFirstOfAlleles()
	{
		ArrayList<String> alleles = new ArrayList<String>(3);
		alleles.add("A");
		alleles.add("C");
		alleles.add("T");
		GeneticVariant testInstance2 = testInstance = ReadOnlyGeneticVariant.createVariant(variantMeta, "rs1", 1, "chr1", null,
				alleles, "T");

		ArrayList<String> allelesExpectedOrder = new ArrayList<String>(3);
		alleles.add("T");
		alleles.add("A");
		alleles.add("C");
		assertEquals(testInstance2.getVariantAlleles().getAllelesAsString(), allelesExpectedOrder);
	}
}
