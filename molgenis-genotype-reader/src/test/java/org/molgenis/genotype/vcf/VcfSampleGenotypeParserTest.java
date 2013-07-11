package org.molgenis.genotype.vcf;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertFalse;
import static org.testng.Assert.assertNotNull;
import static org.testng.Assert.assertTrue;

import java.util.Arrays;
import java.util.List;

import org.molgenis.genotype.vcf.VcfSampleGenotype;
import org.molgenis.genotype.vcf.VcfSampleGenotypeParser;
import org.testng.annotations.Test;

public class VcfSampleGenotypeParserTest
{

	@Test
	public void parseDiploidPased()
	{
		VcfSampleGenotype geno = new VcfSampleGenotypeParser("0|1").parse();
		assertNotNull(geno);

		List<String> alleles = Arrays.asList("A", "T");
		List<String> sampleVariants = geno.getSamleVariants(alleles);
		assertNotNull(sampleVariants);
		assertEquals(sampleVariants.size(), 2);
		assertEquals(sampleVariants.get(0), "A");
		assertEquals(sampleVariants.get(1), "T");

		assertNotNull(geno.getPhasing());
		assertEquals(geno.getPhasing().size(), 1);
		assertTrue(geno.getPhasing().get(0));
	}

	@Test
	public void parseTriploidUnphased()
	{
		VcfSampleGenotype geno = new VcfSampleGenotypeParser("0/2/1").parse();
		assertNotNull(geno);

		List<String> alleles = Arrays.asList("A", "CA", "T");
		List<String> sampleVariants = geno.getSamleVariants(alleles);
		assertNotNull(sampleVariants);
		assertEquals(sampleVariants.size(), 3);
		assertEquals(sampleVariants.get(0), "A");
		assertEquals(sampleVariants.get(1), "T");
		assertEquals(sampleVariants.get(2), "CA");

		assertNotNull(geno.getPhasing());
		assertEquals(geno.getPhasing().size(), 2);
		assertFalse(geno.getPhasing().get(0));
		assertFalse(geno.getPhasing().get(1));
	}
}
