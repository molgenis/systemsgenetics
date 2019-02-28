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
		assertEquals(ZScores.zToP(10), 1.523971e-23, 0.00000001);
		assertEquals(ZScores.zToP(-10), 1.523971e-23, 0.00000001);
		assertEquals(ZScores.zToP(0), 1, 0.00000001);
		assertEquals(ZScores.zToP(2), 0.04550026, 0.00000001);
		assertEquals(ZScores.zToP(-2), 0.04550026, 0.00000001);
	}

	@Test
	public void pToZTwoTailed() {

		assertEquals(ZScores.pToZTwoTailed(0), -40, 0.00000001);
		assertEquals(ZScores.pToZTwoTailed(1), 0, 0.00000001);
		assertEquals(ZScores.pToZTwoTailed(1.523971e-23), -10, 0.00001);
		assertEquals(ZScores.pToZTwoTailed(1e-308), -37.55912122001427, 0.00001);
		assertEquals(ZScores.pToZTwoTailed(0.04550026), -2, 0.00001);
		
	}

}
