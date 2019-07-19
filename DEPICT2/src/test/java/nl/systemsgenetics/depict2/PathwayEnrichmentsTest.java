/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

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
