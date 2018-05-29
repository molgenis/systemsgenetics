package test.java;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import main.java.decon_eQTL.Statistics;

public class StatisticsTest {
	Statistics statistics;
	
	@Before
	public void init() {
	}
	
	@After
	public void tearDown() throws Exception {
		//deleteDir(new File(outputDir));
	}	

	@Test
	public void calculateSpearmanTwoTailedPvalueTest() throws Exception {
		double expectedPval = 1.773084e-06;
		double observedPval = Statistics.calculateSpearmanTwoTailedPvalue(0.1504121, 1000);
		assertEquals(expectedPval,observedPval, 0.001);
		
		expectedPval = 1;
		observedPval = Statistics.calculateSpearmanTwoTailedPvalue(0, 1000);
		assertEquals(expectedPval,observedPval, 0.001);	
	}
	
	@Test
	public void calculateSpearmanTwoTailedPvalueNaNTest() throws Exception {
		double expectedPval = 0;
		double observedPval = Statistics.calculateSpearmanTwoTailedPvalue(1000000, 1000);
		assertEquals(expectedPval,observedPval, 0.001);	
	}
	
	@Test
	public void anovaTest() throws Exception {
		double expectedPval = 4.09E-7;	
		double observedPval = Statistics.anova(31, 22, 100, 101, true);
		assertEquals(expectedPval, observedPval, 0.000000001);
		
		expectedPval = 4.67E-7;	
		observedPval = Statistics.anova(31, 22, 100, 101, false);
		assertEquals(expectedPval, observedPval, 0.000000001);
	}
	
	@Test
	public void anovaErrorTest() throws Exception {
		try {
			Statistics.anova(0, 10, 100, 100, true);
			fail( "My method didn't throw when I expected it to" );
		} catch (RuntimeException expectedException) {
			assertEquals(expectedException.getMessage(), 
						"meanSquareError should not be 0, no variance in the data?");
		}
	}
}