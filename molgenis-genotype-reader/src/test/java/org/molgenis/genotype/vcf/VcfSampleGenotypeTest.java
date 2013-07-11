package org.molgenis.genotype.vcf;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertNotNull;
import static org.testng.Assert.assertNull;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.molgenis.genotype.vcf.VcfSampleGenotype;
import org.testng.annotations.Test;

public class VcfSampleGenotypeTest
{

	@Test
	public void getSamleVariants()
	{
		List<Character> alleleIndices = Arrays.asList('0', '.', '1');
		VcfSampleGenotype sampleGeno = new VcfSampleGenotype(alleleIndices, Collections.<Boolean> emptyList());
		List<String> alleles = Arrays.asList("A", "T");
		List<String> sampleVariants = sampleGeno.getSamleVariants(alleles);
		assertNotNull(sampleVariants);
		assertEquals(sampleVariants.size(), 3);
		assertEquals(sampleVariants.get(0), "A");
		assertNull(sampleVariants.get(1));
		assertEquals(sampleVariants.get(2), "T");
	}

}
