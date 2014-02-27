package org.molgenis.genotype.vcf;

import static org.testng.Assert.assertEquals;

import java.io.IOException;
import java.net.URISyntaxException;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.ResourceTest;
import org.molgenis.genotype.util.CalledDosageConvertor;
import org.molgenis.genotype.variant.GeneticVariant;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

public class VcfGenotypeDataTest2 extends ResourceTest
{
	private VcfGenotypeData genotypeData;

	@BeforeClass
	public void beforeClass() throws IOException, URISyntaxException
	{
		genotypeData = new VcfGenotypeData(getTestVcfGz2(), getTestVcfGz2Tbi());
	}
	
	@Test
	public void getDosage(){
		
		GeneticVariant var = genotypeData.getSnpVariantByPos("21", 9825966);
		
		assertEquals(var.getVariantAlleles(), Alleles.createAlleles(Allele.A, Allele.C));
		
		assertEquals(var.getSampleVariants().get(0), Alleles.createAlleles(Allele.A, Allele.A));
		
		var = genotypeData.getSnpVariantByPos("21", 9834467);
		
		assertEquals(var.getVariantAlleles(), Alleles.createAlleles(Allele.A, Allele.G));
		assertEquals(var.getSampleVariants().get(2), Alleles.createAlleles(Allele.G, Allele.A));
		
		
		
		assertEquals(var.getSampleDosages()[0],1.878, 0.01);
		assertEquals(var.getSampleDosages()[2],1.593, 0.01);
		
	}
	
	@Test
	public void getSamples(){
	
		assertEquals(genotypeData.getSamples().get(2).getId(), "DRS001142");
		
	}


}
