package org.molgenis.genotype.vcf;

import static org.testng.Assert.assertEquals;

import java.io.IOException;
import java.net.URISyntaxException;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.ResourceTest;
import org.molgenis.genotype.variant.GeneticVariant;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

public class VcfGenotypeDataTest2 extends ResourceTest
{
	private VcfGenotypeData genotypeData;

	@BeforeClass
	public void beforeClass() throws IOException, URISyntaxException
	{
		genotypeData = new VcfGenotypeData(getTestVcfGz2(), getTestVcfGz2Tbi(), 0.8);
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

	@Test
	public void getProbs(){
		
		GeneticVariant var = genotypeData.getSnpVariantByPos("21", 9834467);
		
		assertEquals(var.getVariantAlleles(), Alleles.createAlleles(Allele.A, Allele.G));
		assertEquals(var.getSampleVariants().get(2), Alleles.createAlleles(Allele.G, Allele.A));
		
		assertEquals(var.getSampleGenotypeProbilities()[0][0],0.878, 0.01);
		assertEquals(var.getSampleGenotypeProbilities()[0][1],0.122, 0.01);
		assertEquals(var.getSampleGenotypeProbilities()[0][2],0, 0.01);
		
		assertEquals(var.getSampleGenotypeProbilities()[2][0],0.635, 0.01);
		assertEquals(var.getSampleGenotypeProbilities()[2][1],0.323, 0.01);
		assertEquals(var.getSampleGenotypeProbilities()[2][2],0.042, 0.01);
		
	}
	

}
