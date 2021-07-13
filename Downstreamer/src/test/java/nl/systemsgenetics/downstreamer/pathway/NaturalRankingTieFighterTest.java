/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.pathway;

import nl.systemsgenetics.downstreamer.pathway.NaturalRankingTieFighter;
import java.util.Arrays;
import org.apache.commons.math3.stat.ranking.NaNStrategy;
import org.apache.commons.math3.stat.ranking.TiesStrategy;
import static org.testng.Assert.*;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

/**
 *
 * @author patri
 */
public class NaturalRankingTieFighterTest {
	
	public NaturalRankingTieFighterTest() {
	}

	/**
	 * Test of rank method, of class NaturalRankingTieFighter.
	 */
	@Test
	public void testRank() {
		System.out.println("rank");
		double[] data = new double[]{1,2,2,1.5,4,4,6,6,6};
		double[] xWings = new double[]{1,2,3,1,4,4,4,2,6};
		NaturalRankingTieFighter instance = new NaturalRankingTieFighter(NaNStrategy.FAILED,
				TiesStrategy.AVERAGE);
		
		double[] expResult = new double[]{1,3,4,2,5.5,5.5,8,7,9};
		double[] result = instance.rank(data, xWings);
		
		assertTrue(Arrays.equals(expResult, result));
	}
	
}
