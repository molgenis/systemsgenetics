package decon_eQTLTests;

import static org.junit.Assert.assertEquals;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import decon_eQTL.Statistics;

public class StatisticsTest {
	Statistics statistics;
	
	@Before
	public void init() {
		statistics = new Statistics();
	}
	
	@After
	public void tearDown() throws Exception {
		//deleteDir(new File(outputDir));
	}	

}