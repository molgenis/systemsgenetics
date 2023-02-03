/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.runners;

import java.util.ArrayList;
import java.util.List;
import nl.systemsgenetics.downstreamer.DownstreamerStep2Results;
import nl.systemsgenetics.downstreamer.pathway.PathwayDatabase;
import nl.systemsgenetics.downstreamer.runners.options.DownstreamerMode;
import nl.systemsgenetics.downstreamer.runners.options.OptionsModeEnrichment;
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
public class DownstreamerEnrichmentTest {
	
	public DownstreamerEnrichmentTest() {
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
	 * Test of enrichmentAnalysis method, of class DownstreamerEnrichment.
	 */
	@Test
	public void testEnrichmentAnalysis() throws Exception {
		
		
		List<PathwayDatabase> pathwayDatabases = new ArrayList<>();
		pathwayDatabases.add(new PathwayDatabase("test", location, false));
	
		
		OptionsModeEnrichment options = new OptionsModeEnrichment(null, pathwayDatabases, true, true, true, geneInfoFile, singleGwasFile, gwasPvalueMatrixPath, true, true, geneGeneCorrelationPrefix, 0, outputBasePath, logFile, DownstreamerMode.STEP1, true, true);
		DownstreamerStep2Results expResult = null;
		DownstreamerStep2Results result = DownstreamerEnrichment.enrichmentAnalysis(options);
		assertEquals(result, expResult);
		// TODO review the generated test code and remove the default call to fail.
		fail("The test case is a prototype.");
	}
	
}
