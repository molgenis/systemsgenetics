/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats;

import static org.testng.Assert.*;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

/**
 *
 * @author patri
 */
public class ZScoresTest {
	
	public ZScoresTest() {
	}

	@BeforeClass
	public static void setUpClass() throws Exception {
	}

	@BeforeMethod
	public void setUpMethod() throws Exception {
	}

	/**
	 * Test of zToP method, of class ZScores.
	 */
	@Test
	public void testZToP() {
		System.out.println("zToP");
		assertEquals(ZScores.zToP(10), 1.523971e-23, 0.00000001);
		assertEquals(ZScores.zToP(-10), 1.523971e-23, 0.00000001);
		assertEquals(ZScores.zToP(0), 1, 0.00000001);
		assertEquals(ZScores.zToP(2), 0.04550026, 0.00000001);
		assertEquals(ZScores.zToP(-2), 0.04550026, 0.00000001);
	}

	
}
