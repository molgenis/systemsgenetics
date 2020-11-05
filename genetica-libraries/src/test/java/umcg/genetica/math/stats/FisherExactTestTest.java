/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats;

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
public class FisherExactTestTest {
	
	public FisherExactTestTest() {
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
	 * Test of getFisherPValue method, of class FisherExactTest.
	 */
	@Test
	public void testGetFisherPValue() {
		System.out.println("getFisherPValue");
		int n11 = 18030 ;
		int n12 = 1078;
		int n21 = 386;
		int n22 = 22;
		FisherExactTest instance = new FisherExactTest();
		double expResult = 0.9137707;
		double result = instance.getFisherPValue(n11, n12, n21, n22);
		assertEquals(result, expResult, 0.00001);
		
		
		n11 =  19010  ;
		n12 =387 ;
		n21 = 98 ;
		n22 = 21;
		instance = new FisherExactTest();
		expResult = 5.648946e-14;
		result = instance.getFisherPValue(n11, n12, n21, n22);
		assertEquals(result, expResult, 0.0000000000001);
	}

	
}
