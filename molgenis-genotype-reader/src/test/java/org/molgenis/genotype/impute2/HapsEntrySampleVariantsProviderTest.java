package org.molgenis.genotype.impute2;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertFalse;
import static org.testng.Assert.assertTrue;

import java.util.Arrays;
import java.util.List;

import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.ReadOnlyGeneticVariant;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

public class HapsEntrySampleVariantsProviderTest
{
	private HapsEntrySampleVariantsProvider hapsEntrySampleVariantsProvider;
	private GeneticVariant variant;

	@BeforeClass
	public void setUp()
	{
		String line = "7	SNP1	123	A	G	0	0	1*	0*	0	0	1	1";
		HapsEntry entry = HapsLineParser.parse(line);
		hapsEntrySampleVariantsProvider = new HapsEntrySampleVariantsProvider(entry);
		variant = ReadOnlyGeneticVariant.createSnp("SNP1", 123, "7", hapsEntrySampleVariantsProvider, 'A', 'G');
	}

	@Test
	public void getSamplePhasing()
	{
		List<Boolean> phasing = hapsEntrySampleVariantsProvider.getSamplePhasing(variant);
		assertEquals(phasing.size(), 4);
		assertTrue(phasing.get(0));
		assertFalse(phasing.get(1));
		assertTrue(phasing.get(2));
		assertTrue(phasing.get(3));
	}

	@Test
	public void getSampleVariants()
	{
		List<Alleles> samples = hapsEntrySampleVariantsProvider.getSampleVariants(variant);
		assertEquals(
				samples,
				Arrays.asList(Alleles.createBasedOnChars('A', 'A'), Alleles.createBasedOnChars('G', 'A'),
						Alleles.createBasedOnChars('A', 'A'), Alleles.createBasedOnChars('G', 'G')));
	}
}
