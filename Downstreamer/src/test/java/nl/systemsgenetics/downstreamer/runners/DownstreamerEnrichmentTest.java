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
import static nl.systemsgenetics.downstreamer.Downstreamer.DATE_TIME_FORMAT;
import static nl.systemsgenetics.downstreamer.Downstreamer.initializeLoggers;
import nl.systemsgenetics.downstreamer.DownstreamerStep2Results;
import nl.systemsgenetics.downstreamer.gene.Gene;
import nl.systemsgenetics.downstreamer.io.ExcelWriter;
import nl.systemsgenetics.downstreamer.io.IoUtils;
import nl.systemsgenetics.downstreamer.pathway.PathwayDatabase;
import nl.systemsgenetics.downstreamer.runners.options.DownstreamerMode;
import nl.systemsgenetics.downstreamer.runners.options.OptionsModeEnrichment;
import nl.systemsgenetics.downstreamer.summarystatistic.LinearRegressionResult;
import org.apache.logging.log4j.Level;
import static org.testng.Assert.*;
import org.testng.annotations.AfterClass;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import static umcg.genetica.math.matrix2.DoubleMatrixDataset.loadDoubleTextData;

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

		File pathwayFile = new File(this.getClass().getResource("/random/random_pathways.datg").toURI());
		File genesFile = new File(this.getClass().getResource("/random/genesMatchingTest.txt").toURI());
		File gwasFile = new File(this.getClass().getResource("/random/random_gwas_incSnpMinP.tsv").toURI());
		String geneGeneCorrelationPrefix = new File(gwasFile.getParentFile(), "corPerArm").getAbsolutePath() + "/genecor_chr_";
		
		File expectedBetaFile = new File(this.getClass().getResource("/random/random_ds2_sigma_x_df_107_betas").toURI());

//		List<Gene> genes = IoUtils.readGenes(genesFile);
//		ArrayList<String> geneNames = new ArrayList<>();
//		for (Gene gene : genes) {
//			geneNames.add(gene.getGene());
//		}
//		DoubleMatrixDataset<String, String> pathways = DoubleMatrixDataset.loadDoubleBinaryData(pathwayFile.getAbsolutePath());
//		DoubleMatrixDataset<String, String> gwas = DoubleMatrixDataset.loadDoubleData(gwasFile.getAbsolutePath());
//		geneNames.retainAll(pathways.getRowObjects());
//		geneNames.retainAll(gwas.getRowObjects());
//		pathways.viewRowSelection(geneNames).saveBinary(pathwayFile.getParentFile().getAbsolutePath() + "/random_pathwaysNoHla");
//		gwas.viewRowSelection(geneNames).save(pathwayFile.getParentFile().getAbsolutePath() + "/random_gwas_incSnpMinPNoHla.tsv");

		File outdir = getTmpDir();
		
		System.out.println("TMPdir: " + outdir);

		File outputBasePath = new File(outdir.getAbsolutePath() + "/testRun");

		File logFile = new File(outputBasePath.getAbsolutePath() + ".log");

		String startDateTime = DATE_TIME_FORMAT.format(new Date());
		initializeLoggers(Level.INFO, logFile, startDateTime);

		List<PathwayDatabase> pathwayDatabases = new ArrayList<>();
		pathwayDatabases.add(new PathwayDatabase("test", pathwayFile.getAbsolutePath(), false));

		OptionsModeEnrichment options = new OptionsModeEnrichment(null, pathwayDatabases, false, false, false, genesFile, gwasFile, null, false, false, geneGeneCorrelationPrefix, 1, outputBasePath, logFile, DownstreamerMode.STEP1, true, true, true);
		
		options.printOptions();
		
		DownstreamerStep2Results results = DownstreamerEnrichment.enrichmentAnalysis(options);
		
		LinearRegressionResult firstTraitEnrichmentRes = DownstreamerEnrichment.getFirstRestult();
		
		DoubleMatrixDataset<String, String> resBetas = firstTraitEnrichmentRes.getBeta();
		
		DoubleMatrixDataset<String, String> expectedBetas = DoubleMatrixDataset.loadDoubleTextData(expectedBetaFile.getAbsolutePath(), '\t');
		
		DownstreamerRegressionEngineTest.compareTwoMatrices(resBetas, expectedBetas, 0.000000000001);
		
		
		DoubleMatrixDataset<String, String> betas = results.getPathwayEnrichments().get(0).getBetas();
		betas.save(outputBasePath + "_betas.txt");
		
		ExcelWriter writer = new ExcelWriter(results.getGenePvalues().getColObjects(), options);
		writer.saveStep2Excel(results);
	}

	private File getTmpDir() {
		File tmpDir = new File(System.getProperty("java.io.tmpdir"));

		DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
		Date date = new Date();

		//File tmpOutputFolder = new File(tmpDir, "DownstreamerTest" + dateFormat.format(date));
		File tmpOutputFolder = new File(tmpDir, "DownstreamerTestWithHla");

//		Runtime.getRuntime().addShutdownHook(new Thread() {
//			@Override
//			public void run() {
//				if (tmpOutputFolder.isDirectory()) {
//					System.out.println("Removing tmp dir and files");
//					for (File file : tmpOutputFolder.listFiles()) {
//						System.out.println(" - Deleting: " + file.getAbsolutePath());
//						file.delete();
//					}
//					System.out.println(" - Deleting: " + tmpOutputFolder.getAbsolutePath());
//					tmpOutputFolder.delete();
//				}
//			}
//		});
		tmpOutputFolder.mkdir();

		return tmpOutputFolder;
	}

}
