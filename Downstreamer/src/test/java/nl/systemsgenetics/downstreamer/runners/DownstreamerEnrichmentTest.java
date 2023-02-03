/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.runners;

import java.io.File;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
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

		File pathwayFile = new File(this.getClass().getResource("/random/random_pathways.tsv.gz").toURI());
		File genesFile = new File(this.getClass().getResource("/random/genes.txt").toURI());
		File gwasFile = new File(this.getClass().getResource("/random/random_gwas_incSnpMinP.tsv").toURI());
		String geneGeneCorrelationPrefix = new File(gwasFile.getParentFile(), "corPerArm").getAbsolutePath() + "/genecor_chr_";

		File outdir = getTmpDir();
		
		File outputBasePath = new File(outdir.getAbsolutePath() + "testRun");
		
		File logFile = new File(outputBasePath.getAbsolutePath() + ".log");

		List<PathwayDatabase> pathwayDatabases = new ArrayList<>();
		pathwayDatabases.add(new PathwayDatabase("test", pathwayFile.getAbsolutePath(), false));

		OptionsModeEnrichment options = new OptionsModeEnrichment(null, pathwayDatabases, false, false, false, genesFile, gwasFile, null, true, false, geneGeneCorrelationPrefix, 0, outputBasePath, logFile, DownstreamerMode.STEP1, true, true);
		DownstreamerStep2Results expResult = null;
		DownstreamerStep2Results result = DownstreamerEnrichment.enrichmentAnalysis(options);
		assertEquals(result, expResult);
		// TODO review the generated test code and remove the default call to fail.
		fail("The test case is a prototype.");
	}

	private File getTmpDir() {
		File tmpDir = new File(System.getProperty("java.io.tmpdir"));

		DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
		Date date = new Date();

		File tmpOutputFolder = new File(tmpDir, "DoubleMatrixDatasetTest" + dateFormat.format(date));

		Runtime.getRuntime().addShutdownHook(new Thread() {
			@Override
			public void run() {
				if (tmpOutputFolder.isDirectory()) {
					System.out.println("Removing tmp dir and files");
					for (File file : tmpOutputFolder.listFiles()) {
						System.out.println(" - Deleting: " + file.getAbsolutePath());
						file.delete();
					}
					System.out.println(" - Deleting: " + tmpOutputFolder.getAbsolutePath());
					tmpOutputFolder.delete();
				}
			}
		});

		tmpOutputFolder.mkdir();

		return tmpDir;
	}

}
