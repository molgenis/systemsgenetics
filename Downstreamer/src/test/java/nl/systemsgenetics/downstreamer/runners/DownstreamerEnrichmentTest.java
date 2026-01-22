/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.runners;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DoubleStatistic;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import cern.jet.stat.Descriptive;
import java.io.File;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.List;
import static nl.systemsgenetics.downstreamer.Downstreamer.DATE_TIME_FORMAT;
import static nl.systemsgenetics.downstreamer.Downstreamer.initializeLoggers;
import nl.systemsgenetics.downstreamer.DownstreamerStep2Results;
import nl.systemsgenetics.downstreamer.io.ExcelWriter;
import nl.systemsgenetics.downstreamer.pathway.PathwayDatabase;
import nl.systemsgenetics.downstreamer.runners.options.DownstreamerMode;
import nl.systemsgenetics.downstreamer.runners.options.OptionsModeEnrichment;
import nl.systemsgenetics.downstreamer.summarystatistic.LinearRegressionResult;
import org.apache.commons.math3.distribution.TDistribution;
import org.apache.logging.log4j.Level;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

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

		final File pathwayFile = new File(this.getClass().getResource("/random/300G_x_pathways.datg").toURI());
		final File genesFile = new File(this.getClass().getResource("/random/300G_gene_info.txt").toURI());
		final File gwasFile = new File(this.getClass().getResource("/random/300G_y_gwas2.txt").toURI());
		final String geneGeneCorrelationPrefix = new File(gwasFile.getParentFile(), "corPerArm").getAbsolutePath() + "/genecor_chr_";
	
		final File expectedBetaAndSeFile = new File(this.getClass().getResource("/random/300G_beta_se_full_r_calcultated_eigen_decomp_no_intercept.txt").toURI());
		final String outputPrefix = "test1";

		doEnrichmentTest(outputPrefix, pathwayFile, genesFile, gwasFile, geneGeneCorrelationPrefix, expectedBetaAndSeFile, false, true);
	}
	
	@Test
	public void testEnrichmentAnalysis2() throws Exception {

		final File pathwayFile = new File(this.getClass().getResource("/random/300G_x_pathways2.datg").toURI());
		final File genesFile = new File(this.getClass().getResource("/random/300G_gene_info_randomOrder.txt").toURI());
		final File gwasFile = new File(this.getClass().getResource("/random/300G_y_gwas2.txt").toURI());
		final String geneGeneCorrelationPrefix = new File(gwasFile.getParentFile(), "corPerArm").getAbsolutePath() + "/genecor_chr_";
		
		final File expectedBetaAndSeFile = new File(this.getClass().getResource("/random/300G_beta_se_full_r_calcultated_eigen_decomp_no_intercept.txt").toURI());
		final String outputPrefix = "test2";

		doEnrichmentTest(outputPrefix, pathwayFile, genesFile, gwasFile, geneGeneCorrelationPrefix, expectedBetaAndSeFile, false, true);
	}

		@Test
	public void testEnrichmentAnalysis3() throws Exception {

		final File pathwayFile = new File(this.getClass().getResource("/random/300G_x_pathways.datg").toURI());
		final File genesFile = new File(this.getClass().getResource("/random/300G_gene_info.txt").toURI());
		final File gwasFile = new File(this.getClass().getResource("/random/300G_y_gwas2.txt").toURI());
		final String geneGeneCorrelationPrefix = new File(gwasFile.getParentFile(), "corPerArm").getAbsolutePath() + "/genecor_chr_";
		
		final File expectedBetaAndSeFile = new File(this.getClass().getResource("/random/300G_beta_se_full_r_calcultated_eigen_decomp_no_intercept_nochr6.txt").toURI());
		final String outputPrefix = "test3";

		doEnrichmentTest(outputPrefix, pathwayFile, genesFile, gwasFile, geneGeneCorrelationPrefix, expectedBetaAndSeFile, true, true);
	}
	
		@Test
	public void testEnrichmentAnalysis4() throws Exception {

		final File pathwayFile = new File(this.getClass().getResource("/random/300G_x_pathways2.datg").toURI());
		final File genesFile = new File(this.getClass().getResource("/random/300G_gene_info_randomOrder.txt").toURI());
		final File gwasFile = new File(this.getClass().getResource("/random/300G_y_gwas2.txt").toURI());
		final String geneGeneCorrelationPrefix = new File(gwasFile.getParentFile(), "corPerArm").getAbsolutePath() + "/genecor_chr_";
		
		final File expectedBetaAndSeFile = new File(this.getClass().getResource("/random/300G_beta_se_full_r_calcultated_eigen_decomp_no_intercept_nochr6.txt").toURI());
		final String outputPrefix = "test4";

		doEnrichmentTest(outputPrefix, pathwayFile, genesFile, gwasFile, geneGeneCorrelationPrefix, expectedBetaAndSeFile, true, false);
	}
	
	
	private void doEnrichmentTest(final String outputPrefix, final File pathwayFile, final File genesFile, final File gwasFile, final String geneGeneCorrelationPrefix, final File expectedBetaAndSeFile, final boolean excludeHla, final boolean jblas) throws IOException, Exception {
		File outdir = getTmpDir();
		
		System.out.println("TMPdir: " + outdir);

		File outputBasePath = new File(outdir.getAbsolutePath() + "/" + outputPrefix);

		File logFile = new File(outputBasePath.getAbsolutePath() + ".log");

		String startDateTime = DATE_TIME_FORMAT.format(new Date());
		initializeLoggers(Level.INFO, logFile, startDateTime);

		List<PathwayDatabase> pathwayDatabases = new ArrayList<>();
		pathwayDatabases.add(new PathwayDatabase("test", pathwayFile.getAbsolutePath(), false));

		OptionsModeEnrichment options = new OptionsModeEnrichment(null, pathwayDatabases, false, false, false, genesFile, gwasFile, null, excludeHla, true, geneGeneCorrelationPrefix, 1, outputBasePath, logFile, DownstreamerMode.STEP1, true, jblas, true, 0.05);
		
		options.printOptions();
		
		DownstreamerStep2Results results = DownstreamerEnrichment.enrichmentAnalysis(options);
		
		LinearRegressionResult firstTraitEnrichmentRes = DownstreamerEnrichment.getFirstRestult();
		
		DoubleMatrixDataset<String, String> resBetas = firstTraitEnrichmentRes.getBeta();
		DoubleMatrixDataset<String, String> resTs = firstTraitEnrichmentRes.getTstats();
		
		
		DoubleMatrixDataset<String, String> expectedBetasAndSe = DoubleMatrixDataset.loadDoubleTextData(expectedBetaAndSeFile.getAbsolutePath(), '\t');
		

		
		//DownstreamerRegressionEngineTest.compareTwoMatrices(resBetas, expectedBetas, 0.01);
		//resBetas.getCol(0).
		DenseDoubleMatrix2D combinedBetas = new DenseDoubleMatrix2D(resBetas.rows(), 2);
		combinedBetas.viewColumn(0).assign(resBetas.getCol(0));
		combinedBetas.viewColumn(1).assign(expectedBetasAndSe.getCol(0));
		
		DenseDoubleMatrix2D combinedTs = new DenseDoubleMatrix2D(resTs.rows(), 2);
		combinedTs.viewColumn(0).assign(resTs.getCol(0));
		
		for(int r = 0 ; r < combinedTs.rows() ; ++r){
			
			double t = expectedBetasAndSe.getElementQuick(r, 0) / expectedBetasAndSe.getElementQuick(r, 1);
			combinedTs.setQuick(r, 1, t);
			
		}
		
		
		double correlationBetas = DoubleStatistic.correlation(DoubleStatistic.covariance(combinedBetas)).getQuick(0, 1);
		double correlationTs = DoubleStatistic.correlation(DoubleStatistic.covariance(combinedTs)).getQuick(0, 1);
		
		System.out.println("Correlation betas: " + correlationBetas);
		System.out.println("Correlation Ts: " + correlationTs);
		
		//We expect subtle differnce in the beta's and se. Instead of an exact match we check for high correlation
		
		Assert.assertTrue(correlationBetas >= 0.98, "The calculated beta's are not correlated to the expected betas.");
		Assert.assertTrue(correlationTs >= 0.98, "The calculated T's are not correlated to the expected T's.");
				
		DoubleMatrixDataset<String, String> betas = results.getPathwayEnrichments().get(0).getBetas();
		betas.save(outputBasePath + "_betas.txt");
		
		ExcelWriter writer = new ExcelWriter(results.getGenePvalues().getColObjects(), options);
		writer.saveStep2Excel(results);
	}

	private File getTmpDir() {
		File tmpDir = new File(System.getProperty("java.io.tmpdir"));

		DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
		Date date = new Date();

		File tmpOutputFolder = new File(tmpDir, "DownstreamerTest" + dateFormat.format(date));

		Runtime.getRuntime().addShutdownHook(new Thread() {
			@Override
			public void run() {
				if (tmpOutputFolder.isDirectory()) {
					System.out.println("Removing tmp dir and files");
					for (File file : tmpOutputFolder.listFiles()) {
						//System.out.println(" - Deleting: " + file.getAbsolutePath());
						file.delete();
					}
					System.out.println(" - Deleting: " + tmpOutputFolder.getAbsolutePath());
					tmpOutputFolder.delete();
				}
			}
		});
		tmpOutputFolder.mkdir();

		return tmpOutputFolder;
	}

}
