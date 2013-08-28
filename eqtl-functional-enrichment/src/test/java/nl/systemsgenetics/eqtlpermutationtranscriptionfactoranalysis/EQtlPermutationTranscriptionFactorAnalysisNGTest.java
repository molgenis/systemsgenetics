/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlpermutationtranscriptionfactoranalysis;

import java.util.ArrayList;
import org.molgenis.genotype.RandomAccessGenotypeData;
import static org.testng.Assert.*;
import org.testng.annotations.AfterClass;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;
import umcg.genetica.io.trityper.EQTL;

/**
 *
 * @author Matthieu
 */
public class EQtlPermutationTranscriptionFactorAnalysisNGTest {
	
	public EQtlPermutationTranscriptionFactorAnalysisNGTest() {
	}

	@BeforeClass
	public static void setUpClass() throws Exception {
	}

	@AfterClass
	public static void tearDownClass() throws Exception {
	}

	@BeforeMethod
	public void setUpMethod() throws Exception {
	}

	@AfterMethod
	public void tearDownMethod() throws Exception {
	}

	/**
	 * Test of main method, of class EQtlPermutationTranscriptionFactorAnalysis.
	 */
	@Test
	public void testMain() {
		System.out.println("main");
		String[] args = null;
		EQtlPermutationTranscriptionFactorAnalysis.main(args);
		// TODO review the generated test code and remove the default call to fail.
		fail("The test case is a prototype.");
	}

	/**
	 * Test of readEQtlResultData method, of class EQtlPermutationTranscriptionFactorAnalysis.
	 */
	@Test
	public void testReadEQtlResultData() throws Exception {
		System.out.println("readEQtlResultData");
		String eqtlFileLocation = "";
		EQtlPermutationTranscriptionFactorAnalysis instance = new EQtlPermutationTranscriptionFactorAnalysis();
		EQTL[] expResult = null;
		EQTL[] result = instance.readEQtlResultData(eqtlFileLocation);
		assertEquals(result, expResult);
		// TODO review the generated test code and remove the default call to fail.
		fail("The test case is a prototype.");
	}

	/**
	 * Test of readEQtlGenotypeData method, of class EQtlPermutationTranscriptionFactorAnalysis.
	 */
	@Test
	public void testReadEQtlGenotypeData() throws Exception {
		System.out.println("readEQtlGenotypeData");
		String genotypeData = "";
		EQtlPermutationTranscriptionFactorAnalysis instance = new EQtlPermutationTranscriptionFactorAnalysis();
		RandomAccessGenotypeData expResult = null;
		RandomAccessGenotypeData result = instance.readEQtlGenotypeData(genotypeData);
		assertEquals(result, expResult);
		// TODO review the generated test code and remove the default call to fail.
		fail("The test case is a prototype.");
	}

	/**
	 * Test of calculateLd method, of class EQtlPermutationTranscriptionFactorAnalysis.
	 */
	@Test
	public void testCalculateLd() {
		System.out.println("calculateLd");
		//EQtlPermutationTranscriptionFactorAnalysis instance = new EQtlPermutationTranscriptionFactorAnalysis();
		//instance.calculateLd();
		// TODO review the generated test code and remove the default call to fail.
		fail("The test case is a prototype.");
	}

	/**
	 * Test of getEQtlsAsGenotypeSnps method, of class EQtlPermutationTranscriptionFactorAnalysis.
	 */
	/*@Test
	public void testGetEQtlsAsGenotypeSnps() {
		System.out.println("getEQtlsAsGenotypeSnps");
		ArrayList<EQtl> eqtls = null;
		RandomAccessGenotypeData genotypeData = null;
		EQtlPermutationTranscriptionFactorAnalysis instance = new EQtlPermutationTranscriptionFactorAnalysis();
		ArrayList expResult = null;
		ArrayList result = instance.getEQtlsAsGenotypeSnps(eqtls, genotypeData);
		assertEquals(result, expResult);
		// TODO review the generated test code and remove the default call to fail.
		fail("The test case is a prototype.");
	}*/

	/**
	 * Test of performAnalysisStep method, of class EQtlPermutationTranscriptionFactorAnalysis.
	 */
	/*@Test
	public void testPerformAnalysisStep() throws Exception {
		System.out.println("performAnalysisStep");
		ArrayList<EQtl> eqtlData = null;
		int windowSize = 0;
		String outputFilePath = "";
		EQtlPermutationTranscriptionFactorAnalysis instance = new EQtlPermutationTranscriptionFactorAnalysis();
		instance.performAnalysisStep(eqtlData, windowSize, outputFilePath);
		// TODO review the generated test code and remove the default call to fail.
		fail("The test case is a prototype.");
	}*/
}
