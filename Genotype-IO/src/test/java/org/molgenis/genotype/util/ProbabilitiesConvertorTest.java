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
 }