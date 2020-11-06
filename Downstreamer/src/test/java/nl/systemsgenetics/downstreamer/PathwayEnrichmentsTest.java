/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer;

import java.io.File;
import java.net.URISyntaxException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import static org.testng.Assert.*;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.jet.math.tdouble.DoubleFunctions;
import nl.systemsgenetics.downstreamer.pathway.PathwayEnrichments;
import org.apache.commons.math3.distribution.TDistribution;
import org.testng.annotations.Test;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class PathwayEnrichmentsTest {

	private File tmpOutputFolder;

	public PathwayEnrichmentsTest() {
	}

	@org.testng.annotations.BeforeClass
	public static void setUpClass() throws Exception {
	}

	@org.testng.annotations.BeforeMethod
	public void setUpMethod() throws Exception {

		File tmpDir = new File(System.getProperty("java.io.tmpdir"));

		DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
		Date date = new Date();

		tmpOutputFolder = new File(tmpDir, "DoubleMatrixDatasetTest" + dateFormat.format(date));

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

		System.out.println("Temp folder with output of this test: " + tmpOutputFolder.getAbsolutePath());

	}


	/**
	 * Not a proper unit test, more sandbox
	 * @throws Exception
	 */
	@Test
	public void testMergeCorrelationMatrices() throws Exception {
		File invCorFile = new File(this.getClass().getResource("/identity6x6.txt").toURI());
		DoubleMatrixDataset<String, String> geneInvCor = DoubleMatrixDataset.loadDoubleTextData(invCorFile.getAbsolutePath(), '\t');
		invCorFile = new File(this.getClass().getResource("/idenity4x4.txt").toURI());
		DoubleMatrixDataset<String, String> geneInvCor2 = DoubleMatrixDataset.loadDoubleTextData(invCorFile.getAbsolutePath(), '\t');
		invCorFile = new File(this.getClass().getResource("/idenity5x5.txt").toURI());
		DoubleMatrixDataset<String, String> geneInvCor3 = DoubleMatrixDataset.loadDoubleTextData(invCorFile.getAbsolutePath(), '\t');

		List<DoubleMatrixDataset<String, String>> matrices = new ArrayList<>();
		matrices.add(geneInvCor2);
		matrices.add(geneInvCor);
		matrices.add(geneInvCor3);

		DoubleMatrixDataset<String, String> out = PathwayEnrichments.mergeCorrelationMatrices(matrices);
		System.out.println("Debug");


	}


	/**
	 * Not a proper unit test, more sandbox for figuring out the matrix algebra in Java
	 * @throws Exception
	 */
	@Test
	public void testEmpericalPvalues() throws Exception {
		File invCorFile = new File(this.getClass().getResource("/identity6x6.txt").toURI());
		File pathwayFile = new File(this.getClass().getResource("/pathwayGeneScores.txt").toURI());
		File gwasFile = new File(this.getClass().getResource("/gwasGeneScores.txt").toURI());

		DoubleMatrixDataset<String, String> geneZscores = DoubleMatrixDataset.loadDoubleTextData(gwasFile.getAbsolutePath(), '\t');
		DoubleMatrixDataset<String, String> genePathwayZscores = DoubleMatrixDataset.loadDoubleTextData(pathwayFile.getAbsolutePath(), '\t');
		DoubleMatrixDataset<String, String> geneInvCor = DoubleMatrixDataset.loadDoubleTextData(invCorFile.getAbsolutePath(), '\t');

		// Analytical P-values
		// R implementation for reference, can be removed later
		// # Determine SE
		// res       <- y - (x %*% beta)
		// sigma.sqr <- (t(res) %*% Sigi %*% res) / (nrow(x) - ncol(x))
		// se        <- c(sqrt(diag(xtxi))) * c(sqrt(sigma.sqr))

		// # Calculate p
		// tstats <- abs(beta / se)
		// pval <- 2 * pt(tstats, df=nrow(x)-1, lower=F)
		//genePathwayZscores = y;
		//geneZscoresPathwayMatched = x;

		// This implementation is only valid for centered and scaled values, as it does not contain an intercept

		// Determine model residuals

		double beta = 0.1;
		double xtxi = 0.1;
		int n = genePathwayZscores.rows();
		int df = n -1;


		DoubleMatrix1D betaX = geneZscores.getCol(0);
	 	betaX.assign(DoubleFunctions.mult(beta));
		DoubleMatrix1D residuals = genePathwayZscores.getMatrix().viewColumn(0);
		residuals.assign(betaX, DoubleFunctions.minus);

		// TODO: Ugly, but dont know how I can do this better as I need a DoubleMatrix2D for matrix mult
		DoubleMatrix2D residualMatrix = residuals.like2D(n, 1);
		for (int r=0; r < n; ++r) {
			residualMatrix.setQuick(r, 0, residuals.get(r));
		}

		// Determine sigma squared
		DoubleMatrix2D part1 = residualMatrix.like(1, n);
		residualMatrix.zMult(geneInvCor.getMatrix(), part1, 1, 0, true, false);
		DoubleMatrix2D part2 = residualMatrix.like(1, 1);
		part1.zMult(residualMatrix, part2, 1, 0, false, false);
		double sigmaSquared =  part2.get(0,0) / df;

		// Determine Se for beta
		double standardError = Math.sqrt(xtxi) * Math.sqrt(sigmaSquared);

		// Determine pvalue
		double tstatistic = Math.abs(beta / standardError);
		double pvalue = new TDistribution(df).cumulativeProbability(-tstatistic)*2;

		System.out.println("Debug");


	}

	/**
	 * Test of glsFixedInvCor method, of class PathwayEnrichments.
	 */
	@Test
	public void testGlsFixedInvCor() throws URISyntaxException, Exception {
//
//		File identityFile = new File(this.getClass().getResource("/idenity4x4.txt").toURI());
//		File invCorFile = new File(this.getClass().getResource("/invCorMatrix.txt").toURI());
//		File pathwayFile = new File(this.getClass().getResource("/pathwayGeneScores.txt").toURI());
//		File gwasFile = new File(this.getClass().getResource("/gwasGeneScores.txt").toURI());
//
//		DoubleMatrixDataset<String, String> identity = DoubleMatrixDataset.loadDoubleTextData(identityFile.getAbsolutePath(), '\t');
//		DoubleMatrixDataset<String, String> geneZscores = DoubleMatrixDataset.loadDoubleTextData(gwasFile.getAbsolutePath(), '\t');
//		DoubleMatrixDataset<String, String> genePathwayZscores = DoubleMatrixDataset.loadDoubleTextData(pathwayFile.getAbsolutePath(), '\t');
//		DoubleMatrixDataset<String, String> geneInvCor = DoubleMatrixDataset.loadDoubleTextData(invCorFile.getAbsolutePath(), '\t');
//		DoubleMatrixDataset expResult = null;
//
//		identity.printMatrix();
//
//		System.out.println(identity.getMatrix().toString());
//
//		DoubleMatrixDataset result = PathwayEnrichments.glsFixedInvCor(geneZscores, genePathwayZscores, identity);
//
//		result.printMatrix();
//
//		//assertEquals(result, expResult);
	}

	/**
	 * Test of performEnrichmentAnalysis method, of class PathwayEnrichments.
	 */
	@Test
	public void testPerformEnrichmentAnalysis() throws Exception {

//		File identity1qFile = new File(this.getClass().getResource("/identity1q.txt").toURI());
//		File identity1pFile = new File(this.getClass().getResource("/identity1p.txt").toURI());
//		File invCor1qFile = new File(this.getClass().getResource("/invCor1q.txt").toURI());
//		File invCor1pFile = new File(this.getClass().getResource("/invCor1p.txt").toURI());
//		File pathwayFile = new File(this.getClass().getResource("/pathwayGeneScores.txt").toURI());
//		File gwasFile = new File(this.getClass().getResource("/gwasGeneScores.txt").toURI());
//		File gwasFileNull = new File(this.getClass().getResource("/gwasGeneScoresNull.txt").toURI());
//
//		DoubleMatrixDataset<String, String> identity1q = DoubleMatrixDataset.loadDoubleTextData(identity1qFile.getAbsolutePath(), '\t');
//		DoubleMatrixDataset<String, String> identity1p = DoubleMatrixDataset.loadDoubleTextData(identity1pFile.getAbsolutePath(), '\t');
//		DoubleMatrixDataset<String, String> geneZscores = DoubleMatrixDataset.loadDoubleTextData(gwasFile.getAbsolutePath(), '\t');
//		DoubleMatrixDataset<String, String> geneZscoresNullGwas = DoubleMatrixDataset.loadDoubleTextData(gwasFileNull.getAbsolutePath(), '\t');
//		DoubleMatrixDataset<String, String> geneInvCor1q = DoubleMatrixDataset.loadDoubleTextData(invCor1qFile.getAbsolutePath(), '\t');
//		DoubleMatrixDataset<String, String> geneInvCor1p = DoubleMatrixDataset.loadDoubleTextData(invCor1pFile.getAbsolutePath(), '\t');
//
//		Map<String, DoubleMatrixDataset<String, String>> invCorMatrixPerChrArm = new HashMap<>();
//		invCorMatrixPerChrArm.put("1q", identity1q);
//		invCorMatrixPerChrArm.put("1p", identity1p);
//
//		List<PathwayDatabase> pathwayDatabases = new ArrayList<>();
//		pathwayDatabases.add(new PathwayDatabase("test", pathwayFile.getAbsolutePath(), true));
//
//		String outputBasePath = tmpOutputFolder.getAbsolutePath() + "testRun";
//		HashSet<String> hlaGenesToExclude = null;
//
//		//HashMap<PathwayDatabase, DoubleMatrixDataset<String, String>> result = PathwayEnrichments.performEnrichmentAnalysis(geneZscores, geneZscoresNullGwas, invCorMatrixPerChrArm, pathwayDatabases, outputBasePath, hlaGenesToExclude);
//
//		Map<String, DoubleMatrixDataset<String, String>> invCorMatrixPerChrArm2 = new HashMap<>();
//		invCorMatrixPerChrArm2.put("1q", geneInvCor1q);
//		invCorMatrixPerChrArm2.put("1p", geneInvCor1p);
//
		//HashMap<PathwayDatabase, DoubleMatrixDataset<String, String>> result2 = PathwayEnrichments.performEnrichmentAnalysis(geneZscores, geneZscoresNullGwas, invCorMatrixPerChrArm2, pathwayDatabases, outputBasePath, hlaGenesToExclude);

//assertEquals(result, expResult);
		// TODO review the generated test code and remove the default call to fail.
	}

}
