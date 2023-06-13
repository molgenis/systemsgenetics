/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math;

import static org.testng.Assert.*;
import org.testng.annotations.AfterClass;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

/**
 *
 * @author patri
 */
public class JohnsonSuPdfTest {
	
	public JohnsonSuPdfTest() {
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
	 * Test of cumulativeProbability method, of class JohnsonSuPdf.
	 */
	@Test
	public void testCumulativeProbability() {
		double x = -0.5067766;
		double mu = 0.002816414;
		double sigma = 0.129416;
		double nu = -0.01164059;
		double tau = 2.994371;
		double expResult = 0.00035085;
		double result = JohnsonSuPdf.cumulativeProbability(x, mu, sigma, nu, tau);
		assertEquals(result, expResult, 0.0000001);

		
		
		x = -0.9193176;
		mu = -0.0003415547;
		sigma = 0.06248173;
		nu = -0.001747628;
		tau = 3.111419;
		expResult = 9.783203e-13;
		result = JohnsonSuPdf.cumulativeProbability(x, mu, sigma, nu, tau);
		System.out.println(result);
		assertEquals(result, expResult, 1e-12);

		
		
		
		
	}
	
	
}
