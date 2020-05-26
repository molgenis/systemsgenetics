/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.util;

import java.util.ArrayList;
import static org.mockito.Mockito.mock;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.DummySampleVariantsProvider;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.GeneticVariantMeta;
import org.molgenis.genotype.variant.ReadOnlyGeneticVariant;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;
import static org.testng.Assert.*;

/**
 *
 * @author Patrick Deelen
 */
public class GenotypeCountCalculatorTest {
	
	public GenotypeCountCalculatorTest() {
	}

	@BeforeMethod
	public void setUpMethod() throws Exception {
	}

	/**
	 * Test of countGenotypes method, of class GenotypeCountCalculator.
	 */
	@Test
	public void testCountGenotypes() {
		
		GeneticVariant variant;
		ArrayList<Alleles> sampleAlleles;
		SampleVariantsProvider sampleAllelesProvider;
		GeneticVariantMeta variantMeta = mock(GeneticVariantMeta.class);
				
		sampleAlleles = new ArrayList<Alleles>();
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('C', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		sampleAllelesProvider = new DummySampleVariantsProvider(sampleAlleles);	
		variant = ReadOnlyGeneticVariant.createVariant(variantMeta, "rs", 1, "chr", sampleAllelesProvider, Alleles.createBasedOnChars(new char[]{'A', 'C'}));
		
		ArrayList<GenotypeCountCalculator.GenotypeCount> counts = GenotypeCountCalculator.countGenotypes(variant);
		
		assertEquals(counts.get(0).getCount(), 2);
		assertEquals(counts.get(0).getGenotype(), Alleles.createAlleles(Allele.A, Allele.A));
		
		assertEquals(counts.get(1).getCount(), 2);
		assertEquals(counts.get(1).getGenotype(), Alleles.createAlleles(Allele.A, Allele.C));
		
		assertEquals(counts.get(2).getCount(), 1);
		assertEquals(counts.get(2).getGenotype(), Alleles.createAlleles(Allele.C, Allele.C));
		
		
	}
	
	@Test
	public void testCountGenotypes2() {
		
		GeneticVariant variant;
		ArrayList<Alleles> sampleAlleles;
		SampleVariantsProvider sampleAllelesProvider;
		GeneticVariantMeta variantMeta = mock(GeneticVariantMeta.class);
				
		sampleAlleles = new ArrayList<Alleles>();
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('C', '0'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'T'));
		sampleAllelesProvider = new DummySampleVariantsProvider(sampleAlleles);	
		variant = ReadOnlyGeneticVariant.createVariant(variantMeta, "rs", 1, "chr", sampleAllelesProvider, Alleles.createBasedOnChars(new char[]{'A', 'C'}));
		
		ArrayList<GenotypeCountCalculator.GenotypeCount> counts = GenotypeCountCalculator.countGenotypes(variant);
		
		assertEquals(counts.get(0).getCount(), 2);
		assertEquals(counts.get(0).getGenotype(), Alleles.createAlleles(Allele.A, Allele.A));
		
		assertEquals(counts.get(1).getCount(), 2);
		assertEquals(counts.get(1).getGenotype(), Alleles.createAlleles(Allele.A, Allele.C));
		
		assertEquals(counts.get(2).getCount(), 0);
		assertEquals(counts.get(2).getGenotype(), Alleles.createAlleles(Allele.C, Allele.C));
		
		assertEquals(counts.get(3).getCount(), 1);
		assertEquals(counts.get(3).getGenotype(), Alleles.createAlleles(Allele.A, Allele.T));
		
		
	}
	
}