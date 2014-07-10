/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats;

import static org.testng.Assert.*;
import org.testng.annotations.Test;

/**
 *
 * @author Patrick Deelen
 */
public class CorrelationTest {
	
	public CorrelationTest() {
	}

	/**
	 * Test of correlateMeanCenteredData method, of class Correlation.
	 */
	@Test
	public void testCorrelateMeanCenteredData_3args() {
		System.out.println("correlateMeanCenteredData");
		double[] x = new double[]{4.5, -0.5, -1.5, -2.5};
		//double[] y = new double[]{0,  2,  1, -3};
		//double[] x = new double[]{10,5,4,3};
		double[] y = new double[]{5,7,6,2};
		double varX = 3.109126;
		double varY = 2.160247;
		double sdXsdY = varX * varY;
		double expResult = 0.2481458;
		double result = Correlation.correlateMeanCenteredData(x, y, sdXsdY);
		assertEquals(result, expResult, 0.000001);
	}

	
}