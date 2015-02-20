/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.binaryInteraction;

import static org.testng.Assert.*;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

/**
 *
 * @author Patrick Deelen
 */
public class BinaryInteractionMetaAnalysisTest {
	
	public BinaryInteractionMetaAnalysisTest() {
	}

	@BeforeClass
	public static void setUpClass() throws Exception {
	}


	@Test
	public void testweightedZscore1() throws Exception {
		
		double[] zscore = {2.3263479, 0.8416212, 0.5244005};
		int[] samples = {1,1,1};
		double expected = 2.13179059;
		double result = BinaryInteractionMetaAnalysis.weightedZscore(zscore, samples);
		assertEquals(result, expected, 0.0000001);
		
	}
	
	@Test
	public void testweightedZscore2() throws Exception {
		
		double[] zscore = {2.3263479, 0.8416212, 0.5244005};
		int[] samples = {10,1,1};
		double expected = 2.438683939;
		double result = BinaryInteractionMetaAnalysis.weightedZscore(zscore, samples);
		assertEquals(result, expected, 0.0000001);
		
	}
	
	@Test
	public void testweightedZscore3() throws Exception {
		
		double[] zscore = {2.3263479, 0.8416212, 0.5244005};
		int[] samples = {10,1,0};
		double expected = 2.39854709;
		double result = BinaryInteractionMetaAnalysis.weightedZscore(zscore, samples);
		assertEquals(result, expected, 0.0000001);
		
	}
	
	@Test
	public void testweightedZscore4() throws Exception {
		
		double[] zscore = {2.3263479};
		int[] samples = {10};
		double expected = 2.3263479;
		double result = BinaryInteractionMetaAnalysis.weightedZscore(zscore, samples);
		assertEquals(result, expected, 0.0000001);
		
	}
}