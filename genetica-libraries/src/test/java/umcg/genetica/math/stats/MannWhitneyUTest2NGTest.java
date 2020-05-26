/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats;

import static org.testng.Assert.*;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

/**
 *
 * @author patri
 */
public class MannWhitneyUTest2NGTest {

	public MannWhitneyUTest2NGTest() {
	}

	private MannWhitneyUTest2 uTest;

	@BeforeMethod
	public void setUpMethod() throws Exception {
		uTest = new MannWhitneyUTest2();
	}

	@Test
	public void test1() {
		double[] x = {10d, 2d, 4d};
		double[] y = {100d, 3d, 5d, 6d};
		uTest.setData(x, y);
		assertEquals(uTest.getU(), 8d, 0.000001, "Error calculating U value");
		assertEquals(uTest.getZ(), 0.7071067811865475d, 0.000001, "Error calculating Z-score");
		assertEquals(uTest.getP(), 0.4795001221869535, 0.000001, "Error calculating p-value");

	}

	@Test
	public void test2() {
		double[] x = {10d, 2d, 4d};
		double[] y = {100d, 3d, 5d, 6d};
		uTest.setData(y, x);
		assertEquals(uTest.getU(), 8d, 0.000001, "Error calculating U value");
		assertEquals(uTest.getZ(), -0.7071067811865475d, 0.000001, "Error calculating Z-score");
		assertEquals(uTest.getP(), 0.4795001221869535, 0.000001, "Error calculating p-value");

	}

}
