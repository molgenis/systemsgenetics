package org.molgenis.genotype.util;

import static org.testng.AssertJUnit.assertEquals;

import java.util.ArrayList;

import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.testng.annotations.Test;

public class MafCalculatorTest
{

	@Test
	public void calculateMaf()
	{

		ArrayList<Alleles> sampleAlleles;
		Allele refAllele;
		Alleles variantAlleles;

		MafResult result;

		sampleAlleles = new ArrayList<Alleles>();
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'T'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'T'));
		variantAlleles = Alleles.createBasedOnChars('A', 'T');
		result = MafCalculator.calculateMaf(variantAlleles, null, sampleAlleles);
		assertEquals(result.getMinorAllele(), Allele.A);
		assertEquals(result.getFreq(), 0.5d, 0.00000001);

		sampleAlleles = new ArrayList<Alleles>();
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'T'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'T'));
		variantAlleles = Alleles.createBasedOnChars('A', 'T');
		refAllele = Allele.T;
		result = MafCalculator.calculateMaf(variantAlleles, refAllele, sampleAlleles);
		assertEquals(result.getMinorAllele(), Allele.T);
		assertEquals(result.getFreq(), 0.5d, 0.00000001);

		sampleAlleles = new ArrayList<Alleles>();
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'T'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', '0'));
		variantAlleles = Alleles.createBasedOnChars('A', 'T');
		refAllele = Allele.T;
		result = MafCalculator.calculateMaf(variantAlleles, refAllele, sampleAlleles);
		assertEquals(result.getMinorAllele(), Allele.T);
		assertEquals(result.getFreq(), 0.333333333333d, 0.01);

	}
}
