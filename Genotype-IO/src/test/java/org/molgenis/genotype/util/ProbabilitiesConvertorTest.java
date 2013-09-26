/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.util;

import java.util.ArrayList;
import java.util.List;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.probabilities.SampleVariantProbabilities;
import static org.testng.Assert.*;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

/**
 *
 * @author Patrick Deelen
 */
public class ProbabilitiesConvertorTest {
	
	public ProbabilitiesConvertorTest() {
	}

	@BeforeMethod
	public void setUpMethod() throws Exception {
	}

	/**
	 * Test of convertCalledAllelesToDosage method, of class ProbabilitiesConvertor.
	 */
	@Test
	public void testConvertCalledAllelesToDosage1() {
		System.out.println("convertCalledAllelesToDosage");
		List<Alleles> sampleAlleles = new ArrayList<Alleles>();
		sampleAlleles.add(Alleles.createAlleles(Allele.A, Allele.A));
		sampleAlleles.add(Alleles.createAlleles(Allele.A, Allele.A));
		sampleAlleles.add(Alleles.createAlleles(Allele.A, Allele.C));
		sampleAlleles.add(Alleles.createAlleles(Allele.C, Allele.A));
		sampleAlleles.add(Alleles.createAlleles(Allele.C, Allele.C));
		sampleAlleles.add(Alleles.createAlleles(Allele.A, Allele.C));
		
		Alleles alleles = Alleles.createAlleles(Allele.A, Allele.C);
		SampleVariantProbabilities[] expResult = new SampleVariantProbabilities[6];
		
		expResult[0] = SampleVariantProbabilities.AA_PROB;
		expResult[1] = SampleVariantProbabilities.AA_PROB;
		expResult[2] = SampleVariantProbabilities.AB_PROB;
		expResult[3] = SampleVariantProbabilities.AB_PROB;
		expResult[4] = SampleVariantProbabilities.BB_PROB;
		expResult[5] = SampleVariantProbabilities.AB_PROB;
		
		SampleVariantProbabilities[] result = ProbabilitiesConvertor.convertCalledAllelesToDosage(sampleAlleles, alleles);
		assertEqualsSampleVariantProbabilitiesArray(result, expResult);
		
	}
	
	@Test
	public void testConvertCalledAllelesToDosage2() {
		System.out.println("convertCalledAllelesToDosage");
		List<Alleles> sampleAlleles = new ArrayList<Alleles>();
		sampleAlleles.add(Alleles.createAlleles(Allele.A, Allele.A));
		sampleAlleles.add(Alleles.createAlleles(Allele.ZERO, Allele.A));
		sampleAlleles.add(Alleles.createAlleles(Allele.A, Allele.C));
		sampleAlleles.add(Alleles.createAlleles(Allele.C, Allele.G));
		sampleAlleles.add(Alleles.createAlleles(Allele.C, Allele.C));
		sampleAlleles.add(Alleles.createAlleles(Allele.A, Allele.C));
		
		Alleles alleles = Alleles.createAlleles(Allele.C, Allele.A);
		SampleVariantProbabilities[] expResult = new SampleVariantProbabilities[6];
		
		expResult[0] = SampleVariantProbabilities.BB_PROB;
		expResult[1] = SampleVariantProbabilities.MISSING_PROB;
		expResult[2] = SampleVariantProbabilities.AB_PROB;
		expResult[3] = SampleVariantProbabilities.MISSING_PROB;
		expResult[4] = SampleVariantProbabilities.AA_PROB;
		expResult[5] = SampleVariantProbabilities.AB_PROB;
		
		SampleVariantProbabilities[] result = ProbabilitiesConvertor.convertCalledAllelesToDosage(sampleAlleles, alleles);
		assertEqualsSampleVariantProbabilitiesArray(result, expResult);
		
	}
	
	@Test
	public void testConvertCalledAllelesToDosage3() {
		System.out.println("convertCalledAllelesToDosage");
		List<Alleles> sampleAlleles = new ArrayList<Alleles>();
		sampleAlleles.add(Alleles.createAlleles(Allele.A, Allele.A));
		sampleAlleles.add(Alleles.createAlleles(Allele.ZERO, Allele.A));
		sampleAlleles.add(Alleles.createAlleles(Allele.A, Allele.C));
		sampleAlleles.add(Alleles.createAlleles(Allele.C, Allele.G));
		sampleAlleles.add(Alleles.createAlleles(Allele.C, Allele.C));
		sampleAlleles.add(Alleles.createAlleles(Allele.A, Allele.C));
		
		Alleles alleles = Alleles.createAlleles(Allele.C);
		SampleVariantProbabilities[] expResult = new SampleVariantProbabilities[6];
		
		expResult[0] = SampleVariantProbabilities.MISSING_PROB;
		expResult[1] = SampleVariantProbabilities.MISSING_PROB;
		expResult[2] = SampleVariantProbabilities.MISSING_PROB;
		expResult[3] = SampleVariantProbabilities.MISSING_PROB;
		expResult[4] = SampleVariantProbabilities.AA_PROB;
		expResult[5] = SampleVariantProbabilities.MISSING_PROB;
		
		SampleVariantProbabilities[] result = ProbabilitiesConvertor.convertCalledAllelesToDosage(sampleAlleles, alleles);
		assertEqualsSampleVariantProbabilitiesArray(result, expResult);
		
	}
	
	@Test(expectedExceptions = GenotypeDataException.class)  
	public void testConvertCalledAllelesToDosage4() {
		System.out.println("convertCalledAllelesToDosage");
		List<Alleles> sampleAlleles = new ArrayList<Alleles>();
		sampleAlleles.add(Alleles.createAlleles(Allele.A, Allele.A));
		sampleAlleles.add(Alleles.createAlleles(Allele.ZERO, Allele.A));
		sampleAlleles.add(Alleles.createAlleles(Allele.A, Allele.C));
		sampleAlleles.add(Alleles.createAlleles(Allele.C, Allele.G));
		sampleAlleles.add(Alleles.createAlleles(Allele.C, Allele.C));
		sampleAlleles.add(Alleles.createAlleles(Allele.A, Allele.C));
		
		Alleles alleles = Alleles.createAlleles();

		SampleVariantProbabilities[] result = ProbabilitiesConvertor.convertCalledAllelesToDosage(sampleAlleles, alleles);
		
	}
	
	@Test
	public void testConvertCalledAllelesToDosage5() {
		System.out.println("convertCalledAllelesToDosage");
		List<Alleles> sampleAlleles = new ArrayList<Alleles>();
		sampleAlleles.add(Alleles.createAlleles(Allele.A, Allele.A));
		sampleAlleles.add(Alleles.createAlleles(Allele.ZERO, Allele.A));
		sampleAlleles.add(Alleles.createAlleles(Allele.A, Allele.C));
		sampleAlleles.add(Alleles.createAlleles(Allele.C, Allele.G));
		sampleAlleles.add(Alleles.createAlleles(Allele.C, Allele.C));
		sampleAlleles.add(Alleles.createAlleles(Allele.A, Allele.C));
		
		Alleles alleles = Alleles.createAlleles(Allele.C, Allele.A, Allele.G);
		SampleVariantProbabilities[] expResult = new SampleVariantProbabilities[6];
		
		expResult[0] = SampleVariantProbabilities.MISSING_PROB;
		expResult[1] = SampleVariantProbabilities.MISSING_PROB;
		expResult[2] = SampleVariantProbabilities.MISSING_PROB;
		expResult[3] = SampleVariantProbabilities.MISSING_PROB;
		expResult[4] = SampleVariantProbabilities.MISSING_PROB;
		expResult[5] = SampleVariantProbabilities.MISSING_PROB;
		
		SampleVariantProbabilities[] result = ProbabilitiesConvertor.convertCalledAllelesToDosage(sampleAlleles, alleles);
		assertEqualsSampleVariantProbabilitiesArray(result, expResult);
		
	}
	
	
	public static void assertEqualsSampleVariantProbabilitiesArray(SampleVariantProbabilities[] result, SampleVariantProbabilities[] expResult){
		
		assertEquals(result.length, expResult.length, "Different number of probs");
		
		for(int i = 0 ; i < result.length ; ++i){
			assertEqualsSampleVariantProbabilities(result[i], expResult[i], "sample " + i);
		}
		
	}
	
	public static void assertEqualsSampleVariantProbabilities(SampleVariantProbabilities result, SampleVariantProbabilities expResult, String comparisonName){
		float[] resultVector = result.getProbilities();
		float[] expResultVector = expResult.getProbilities();
		
		for(int i = 0 ; i < resultVector.length ; ++i){
			assertEquals(resultVector[i], expResultVector[i], 0.0001f, "Different prob value for element " + i + " in prob vector of " + comparisonName);
		}
		
	}
	
}