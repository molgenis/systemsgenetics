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
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;
import static org.molgenis.genotype.util.AssertExtended.*;

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
	 * Test of convertDosageToProbabilityHeuristic method, of class
	 * ProbabilitiesConvertor.
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
		float[][] expResult = new float[6][3];

		expResult[0] = new float[]{1, 0, 0};
		expResult[1] = new float[]{1, 0, 0};
		expResult[2] = new float[]{0, 1, 0};
		expResult[3] = new float[]{0, 1, 0};
		expResult[4] = new float[]{0, 0, 1};
		expResult[5] = new float[]{0, 1, 0};

		float[][] result = ProbabilitiesConvertor.convertCalledAllelesToProbability(sampleAlleles, alleles);
		assertEquals(result, expResult, 0.001f, "Probs not identical");

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

		float[][] expResult = new float[6][3];

		expResult[0] = new float[]{0, 0, 1};
		expResult[1] = new float[]{0, 0, 0};
		expResult[2] = new float[]{0, 1, 0};
		expResult[3] = new float[]{0, 0, 0};
		expResult[4] = new float[]{1, 0, 0};
		expResult[5] = new float[]{0, 1, 0};


		float[][] result = ProbabilitiesConvertor.convertCalledAllelesToProbability(sampleAlleles, alleles);
		assertEquals(result, expResult, 0.001f, "Probs not identical");

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

		float[][] expResult = new float[6][3];

		expResult[0] = new float[]{0, 0, 0};
		expResult[1] = new float[]{0, 0, 0};
		expResult[2] = new float[]{0, 0, 0};
		expResult[3] = new float[]{0, 0, 0};
		expResult[4] = new float[]{1, 0, 0};
		expResult[5] = new float[]{0, 0, 0};


		float[][] result = ProbabilitiesConvertor.convertCalledAllelesToProbability(sampleAlleles, alleles);
		assertEquals(result, expResult, 0.001f, "Probs not identical");

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

		float[][] expResult = new float[6][3];

		expResult[0] = new float[]{0, 0, 1};
		expResult[1] = new float[]{0, 0, 0};
		expResult[2] = new float[]{0, 1, 0};
		expResult[3] = new float[]{0, 0, 0};
		expResult[4] = new float[]{1, 0, 0};
		expResult[5] = new float[]{0, 1, 0};


		float[][] result = ProbabilitiesConvertor.convertCalledAllelesToProbability(sampleAlleles, alleles);

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

		float[][] expResult = new float[6][3];

		expResult[0] = new float[]{0, 0, 0};
		expResult[1] = new float[]{0, 0, 0};
		expResult[2] = new float[]{0, 0, 0};
		expResult[3] = new float[]{0, 0, 0};
		expResult[4] = new float[]{0, 0, 0};
		expResult[5] = new float[]{0, 0, 0};


		float[][] result = ProbabilitiesConvertor.convertCalledAllelesToProbability(sampleAlleles, alleles);
		assertEquals(result, expResult, 0.001f, "Probs not identical");

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

		float[][] expResult = new float[6][3];

		expResult[0] = new float[]{1, 0, 0};
		expResult[1] = new float[]{1, 0, 0};
		expResult[2] = new float[]{0, 1, 0};
		expResult[3] = new float[]{0, 1, 0};
		expResult[4] = new float[]{0, 0, 1};
		expResult[5] = new float[]{0, 1, 0};


		float[][] result = ProbabilitiesConvertor.convertDosageToProbabilityHeuristic(sampleDosage);
		assertEquals(result, expResult, 0.001f, "Probs not identical");


	}

	@Test
	public void testConvertDosageToProbabilityHeuristic2() {
		System.out.println("convertDosageToProbabilityHeuristic");

		float[] sampleDosage = new float[]{0, 1, 2, -1, 3, 0.5f, 1.5f, 1.75f};

		float[][] expResult = new float[8][3];

		expResult[0] = new float[]{0, 0, 1};
		expResult[1] = new float[]{0, 1, 0};
		expResult[2] = new float[]{1, 0, 0};
		expResult[3] = new float[]{0, 0, 0};
		expResult[4] = new float[]{0, 0, 0};
		expResult[5] = new float[]{0, 0.5f, 0.5f};
		expResult[6] = new float[]{0.5f, 0.5f, 0};
		expResult[7] = new float[]{0.75f, 0.25f, 0};

		float[][] result = ProbabilitiesConvertor.convertDosageToProbabilityHeuristic(sampleDosage);
		assertEquals(result, expResult, 0.001f, "Probs not identical");

	}
	
	@Test
	public void testConvertProbabilitiesToAlleles() {
		System.out.println("convertProbabilitiesToAlleles");
		
		float[][] probs = new float[8][3];

		probs[0] = new float[]{1, 0, 0};
		probs[1] = new float[]{1, 0, 0};
		probs[2] = new float[]{0, 1, 0};
		probs[3] = new float[]{0, 1, 0};
		probs[4] = new float[]{0, 0, 1};
		probs[5] = new float[]{0, 1, 0};
		probs[6] = new float[]{0.3f, 0.3f, 0.3f};
		probs[7] = new float[]{0.3f, 0.3f, 0.4f};
		
		List<Alleles> expectedSampleAlleles = new ArrayList<Alleles>();
		expectedSampleAlleles.add(Alleles.createAlleles(Allele.A, Allele.A));
		expectedSampleAlleles.add(Alleles.createAlleles(Allele.A, Allele.A));
		expectedSampleAlleles.add(Alleles.createAlleles(Allele.A, Allele.C));
		expectedSampleAlleles.add(Alleles.createAlleles(Allele.A, Allele.C));
		expectedSampleAlleles.add(Alleles.createAlleles(Allele.C, Allele.C));
		expectedSampleAlleles.add(Alleles.createAlleles(Allele.A, Allele.C));
		expectedSampleAlleles.add(Alleles.createAlleles(Allele.ZERO, Allele.ZERO));
		expectedSampleAlleles.add(Alleles.createAlleles(Allele.C, Allele.C));
		
		Alleles alleles = Alleles.createAlleles(Allele.A, Allele.C);
		
		List<Alleles> sampleAlleles = ProbabilitiesConvertor.convertProbabilitiesToAlleles(probs, alleles, 0.4);
		
		for(int i = 0 ; i < expectedSampleAlleles.size() ; ++i){
			assertEquals(sampleAlleles.get(i), expectedSampleAlleles.get(i), "sample index " + i);
		}
		
	}
	
	@Test
	public void testConvertProbabilitiesToDosage() {
		System.out.println("convertProbabilitiesToAlleles");
		
		float[][] probs = new float[8][3];

		probs[0] = new float[]{1, 0, 0};
		probs[1] = new float[]{1, 0, 0};
		probs[2] = new float[]{0, 1, 0};
		probs[3] = new float[]{0, 1, 0};
		probs[4] = new float[]{0, 0, 1};
		probs[5] = new float[]{0, 1, 0};
		probs[6] = new float[]{0.3f, 0.3f, 0.3f};
		probs[7] = new float[]{0.3f, 0.3f, 0.4f};
		
		float[] expectedDosage = new float[]{2,2,1,1,0,1,-1,0.9f};
		float[] dosage = ProbabilitiesConvertor.convertProbabilitiesToDosage(probs, 0.4);
		
		assertEquals(dosage, expectedDosage, 0.00001f, "");
		
		
	}

	@Test
	public void testBiallelicBgenProbabilitiesToProbabilities() {
		System.out.println("convertBiallelicBgenProbabilitiesToProbabilities");

		double[][] probs = new double[6][];

		probs[0] = new double[]{1, 0};
		probs[1] = new double[]{1, 0, 0};
		probs[2] = new double[]{0, 1, 0};
		probs[3] = new double[]{0, 1, 0};
		probs[4] = new double[]{0, 0, 1};
		probs[5] = new double[]{1, 0, 0, 0};

		float[][] expectedProbs = new float[6][];

		expectedProbs[0] = new float[]{0, 0, 0};
		expectedProbs[1] = new float[]{1, 0, 0};
		expectedProbs[2] = new float[]{0, 1, 0};
		expectedProbs[3] = new float[]{0, 1, 0};
		expectedProbs[4] = new float[]{0, 0, 1};
		expectedProbs[5] = new float[]{0, 0, 0};

		float[][] actualProbabilities = ProbabilitiesConvertor
				.convertBiallelicBgenProbabilitiesToProbabilities(probs);
		for (int i = 0; i < actualProbabilities.length; i++) {
			float[] actualProbs = actualProbabilities[i];
			assertEquals(actualProbs, expectedProbs[i]);
		}
	}

	@Test
	public void testProbabilitiesToBgenProbabilities() {
		System.out.println("convertProbabilitiesToBgenProbabilities");

		float[][] probs = new float[3][];

		probs[0] = new float[]{0, 0, 1};
		probs[1] = new float[]{1, 0, 0};
		probs[2] = new float[]{0, 1, 0};

		double[][] expectedProbs = new double[3][];

		expectedProbs[0] = new double[]{0 ,0, 1};
		expectedProbs[1] = new double[]{1, 0, 0};
		expectedProbs[2] = new double[]{0, 1, 0};

		double[][] actualProbabilities = ProbabilitiesConvertor
				.convertProbabilitiesToBgenProbabilities(probs);
		for (int i = 0; i < actualProbabilities.length; i++) {
			double[] actualProbs = actualProbabilities[i];
			assertEquals(actualProbs, expectedProbs[i]);
		}
	}
	
 }