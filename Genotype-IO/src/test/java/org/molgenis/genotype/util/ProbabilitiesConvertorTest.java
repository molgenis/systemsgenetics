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
import org.molgenis.genotype.probabilities.SampleVariantProbabilities3Probs;
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
	 * Test of convertDosageToProbabilityHeuristic method, of class ProbabilitiesConvertor.
	 */
	@Test
	public void testConvertCalledAllelesToProbability1() {
		System.out.println("convertCalledAllelesToProbability");
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
		
		SampleVariantProbabilities[] result = ProbabilitiesConvertor.convertCalledAllelesToProbability(sampleAlleles, alleles);
		assertEqualsSampleVariantProbabilitiesArray(result, expResult);
		
	}
	
	@Test
	public void testConvertCalledAllelesToProbability2() {
		System.out.println("convertCalledAllelesToProbability");
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
		
		SampleVariantProbabilities[] result = ProbabilitiesConvertor.convertCalledAllelesToProbability(sampleAlleles, alleles);
		assertEqualsSampleVariantProbabilitiesArray(result, expResult);
		
	}
	
	@Test
	public void testConvertCalledAllelesToProbability3() {
		System.out.println("convertCalledAllelesToProbability");
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
		
		SampleVariantProbabilities[] result = ProbabilitiesConvertor.convertCalledAllelesToProbability(sampleAlleles, alleles);
		assertEqualsSampleVariantProbabilitiesArray(result, expResult);
		
	}
	
	@Test(expectedExceptions = GenotypeDataException.class)  
	public void testConvertCalledAllelesToProbability4() {
		System.out.println("convertCalledAllelesToProbability");
		List<Alleles> sampleAlleles = new ArrayList<Alleles>();
		sampleAlleles.add(Alleles.createAlleles(Allele.A, Allele.A));
		sampleAlleles.add(Alleles.createAlleles(Allele.ZERO, Allele.A));
		sampleAlleles.add(Alleles.createAlleles(Allele.A, Allele.C));
		sampleAlleles.add(Alleles.createAlleles(Allele.C, Allele.G));
		sampleAlleles.add(Alleles.createAlleles(Allele.C, Allele.C));
		sampleAlleles.add(Alleles.createAlleles(Allele.A, Allele.C));
		
		Alleles alleles = Alleles.createAlleles();

		SampleVariantProbabilities[] result = ProbabilitiesConvertor.convertCalledAllelesToProbability(sampleAlleles, alleles);
		
	}
	
	@Test
	public void testConvertCalledAllelesToProbability5() {
		System.out.println("convertCalledAllelesToProbability");
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
		
		SampleVariantProbabilities[] result = ProbabilitiesConvertor.convertCalledAllelesToProbability(sampleAlleles, alleles);
		assertEqualsSampleVariantProbabilitiesArray(result, expResult);
		
	}
	
	@Test
	public void testConvertDosageToProbabilityHeuristic1() {
		System.out.println("convertDosageToProbabilityHeuristic");
		List<Alleles> sampleAlleles = new ArrayList<Alleles>();
		sampleAlleles.add(Alleles.createAlleles(Allele.A, Allele.A));
		sampleAlleles.add(Alleles.createAlleles(Allele.A, Allele.A));
		sampleAlleles.add(Alleles.createAlleles(Allele.A, Allele.C));
		sampleAlleles.add(Alleles.createAlleles(Allele.C, Allele.A));
		sampleAlleles.add(Alleles.createAlleles(Allele.C, Allele.C));
		sampleAlleles.add(Alleles.createAlleles(Allele.A, Allele.C));
		
		Alleles alleles = Alleles.createAlleles(Allele.A, Allele.C);
		
		float[] sampleDosage = CalledDosageConvertor.convertCalledAllelesToDosage(sampleAlleles, alleles, null);
	
		SampleVariantProbabilities[] expResult = new SampleVariantProbabilities[6];
		
		expResult[0] = SampleVariantProbabilities.AA_PROB;
		expResult[1] = SampleVariantProbabilities.AA_PROB;
		expResult[2] = SampleVariantProbabilities.AB_PROB;
		expResult[3] = SampleVariantProbabilities.AB_PROB;
		expResult[4] = SampleVariantProbabilities.BB_PROB;
		expResult[5] = SampleVariantProbabilities.AB_PROB;
		
		SampleVariantProbabilities[] result = ProbabilitiesConvertor.convertDosageToProbabilityHeuristic(sampleDosage);
		assertEqualsSampleVariantProbabilitiesArray(result, expResult);

		
	}
	
	@Test
	public void testConvertDosageToProbabilityHeuristic2() {
		System.out.println("convertDosageToProbabilityHeuristic");
		
		
		float[] sampleDosage = new float[] {0, 1, 2, -1, 3, 0.5f, 1.5f, 1.75f};
		
		SampleVariantProbabilities[] expResult = new SampleVariantProbabilities[8];
		
		expResult[0] = SampleVariantProbabilities.BB_PROB;
		expResult[1] = SampleVariantProbabilities.AB_PROB;
		expResult[2] = SampleVariantProbabilities.AA_PROB;
		expResult[3] = SampleVariantProbabilities.MISSING_PROB;
		expResult[4] = SampleVariantProbabilities.MISSING_PROB;
		expResult[5] = new SampleVariantProbabilities3Probs(new float[] {0,0.5f,0.5f});
		expResult[6] = new SampleVariantProbabilities3Probs(new float[] {0.5f,0.5f, 0});
		expResult[7] = new SampleVariantProbabilities3Probs(new float[] {0.75f,0.25f, 0});
		
		SampleVariantProbabilities[] result = ProbabilitiesConvertor.convertDosageToProbabilityHeuristic(sampleDosage);
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