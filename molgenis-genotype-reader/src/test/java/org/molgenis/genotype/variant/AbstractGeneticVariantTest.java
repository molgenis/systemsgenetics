package org.molgenis.genotype.variant;

import static org.testng.AssertJUnit.assertEquals;

import org.molgenis.genotype.DummySampleVariantsProvider;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

public class AbstractGeneticVariantTest
{

	GeneticVariant variant1;
	GeneticVariant variant2;
	GeneticVariant variant3;
	GeneticVariant variant4;
	GeneticVariant variant5;
	GeneticVariant variant6;
	GeneticVariant variant7;
	GeneticVariant variant8;
	GeneticVariant variant9;
	GeneticVariant variant10;

	@BeforeMethod
	public void setup()
	{

		SampleVariantsProvider provider1 = new DummySampleVariantsProvider(null);
		SampleVariantsProvider provider2 = new DummySampleVariantsProvider(null);

		variant1 = ReadOnlyGeneticVariant.createSnp("rs1", 1, "1", provider1, 'A', 'T');
		variant2 = ReadOnlyGeneticVariant.createSnp("rs2", 1, "1", provider1, 'A', 'T');
		variant3 = ReadOnlyGeneticVariant.createSnp("rs3", 2, "1", provider1, 'A', 'T');
		variant4 = ReadOnlyGeneticVariant.createSnp("rs4", 1, "2", provider1, 'A', 'T');
		variant5 = ReadOnlyGeneticVariant.createSnp("rs5", 1, "1", provider1, 'G', 'T');
		variant6 = ReadOnlyGeneticVariant.createVariant("rs6", 1, "1", provider1, "G", "AA");
		variant7 = ReadOnlyGeneticVariant.createVariant("rs7", 1, "1", provider1, "GG", "T");
		variant8 = ReadOnlyGeneticVariant.createVariant("rs8", 1, "X", provider1, "GG", "T");
		variant9 = ReadOnlyGeneticVariant.createSnp("rs1", 1, "1", provider2, 'A', 'T');
		variant10 = ReadOnlyGeneticVariant.createSnp("rs1", 3, "1", provider2, 'T', 'T');

	}

	@Test
	public void compareTo()
	{
		assertEquals(variant1.compareTo(variant3) < 0, true);
		assertEquals(variant3.compareTo(variant1) > 0, true);
		assertEquals(variant1.compareTo(variant4) < 0, true);
		assertEquals(variant4.compareTo(variant1) > 0, true);
		assertEquals(variant3.compareTo(variant4) < 0, true);
		assertEquals(variant4.compareTo(variant3) > 0, true);
		assertEquals(variant3.compareTo(variant8) < 0, true);
		assertEquals(variant8.compareTo(variant3) > 0, true);
		assertEquals(variant1.compareTo(variant10) < 0, true);
		assertEquals(variant10.compareTo(variant1) > 0, true);
	}

	@Test
	public void equals()
	{

		assertEquals(variant1.equals(variant1), true);
		assertEquals(variant1.equals(variant2), true);
		assertEquals(variant1.equals(variant3), false);
		assertEquals(variant1.equals(variant4), false);
		assertEquals(variant1.equals(variant5), false);
		assertEquals(variant1.equals(variant6), false);
		assertEquals(variant1.equals(variant7), false);
		assertEquals(variant1.equals(variant8), false);
		assertEquals(variant1.equals(variant9), false);
		assertEquals(variant2.equals(variant3), false);
		assertEquals(variant5.equals(variant6), false);
		assertEquals(variant7.equals(variant8), false);

	}
}
