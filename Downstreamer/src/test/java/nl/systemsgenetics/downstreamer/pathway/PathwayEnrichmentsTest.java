/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.pathway;

import java.io.File;
import java.net.URISyntaxException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import nl.systemsgenetics.downstreamer.gene.Gene;
import nl.systemsgenetics.downstreamer.pathway.PathwayEnrichments.MetaGene;
import org.testng.Assert;
import static org.testng.Assert.assertEquals;
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
public class PathwayEnrichmentsTest {

	private static File debugFolder;

	public PathwayEnrichmentsTest() {
	}

	@BeforeClass
	public static void setUpClass() throws Exception {

		File tmpDir = new File(System.getProperty("java.io.tmpdir"));

		DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
		Date date = new Date();

		debugFolder = new File(tmpDir, "DoubleMatrixDatasetTest" + dateFormat.format(date));

		debugFolder.mkdir();

		System.out.println("Temp folder with output of this test: " + debugFolder.getAbsolutePath());

	}

	@AfterClass
	public static void tearDownClass() throws Exception {
		if (debugFolder.isDirectory()) {
			System.out.println("Removing tmp dir and files");
			for (File file : debugFolder.listFiles()) {
				System.out.println(" - Deleting: " + file.getAbsolutePath());
				file.delete();
			}
			System.out.println(" - Deleting: " + debugFolder.getAbsolutePath());
			debugFolder.delete();
		}
	}

	@BeforeMethod
	public void setUpMethod() throws Exception {
	}

	@AfterMethod
	public void tearDownMethod() throws Exception {
	}

	/**
	 * Test of createColumnForceNormalDuplicate method, of class
	 * PathwayEnrichments.
	 */
	@Test
//	public void testGroupCorrelatedGenesPerChrArmAndCollapseDatasetToMetaGenes() throws URISyntaxException, Exception {
//
//		File geneZscoresNullGwasFile = new File(this.getClass().getResource("/geneZscoresNullGwas.txt").toURI());
//		File mergedMetaGenesMeanRefFile = new File(this.getClass().getResource("/mergedMetaGenes1.txt").toURI());
//		File mergedMetaGenesZscoreSumRefFile = new File(this.getClass().getResource("/mergedMetaGenes2.txt").toURI());
//
//		DoubleMatrixDataset<String, String> geneZscoresNullGwas = DoubleMatrixDataset.loadDoubleTextData(geneZscoresNullGwasFile.getPath(), '\t');
//		DoubleMatrixDataset<String, String> mergedMetaGenesMeanRef = DoubleMatrixDataset.loadDoubleTextData(mergedMetaGenesMeanRefFile.getPath(), '\t');
//		DoubleMatrixDataset<String, String> mergedMetaGenesZscoreSumRef = DoubleMatrixDataset.loadDoubleTextData(mergedMetaGenesZscoreSumRefFile.getPath(), '\t');
//
//
//		ArrayList<Gene> genes = new ArrayList();
//
//		genes.add(new Gene("1", "1", 0, 0, "q11"));
//		genes.add(new Gene("2", "1", 0, 0, "q12"));
//		genes.add(new Gene("3", "1", 0, 0, "q22"));
//		genes.add(new Gene("4", "1", 0, 0, "q12"));
//		genes.add(new Gene("5", "1", 0, 0, "p12"));
//		genes.add(new Gene("6", "1", 0, 0, "p12"));
//		genes.add(new Gene("7", "1", 0, 0, "p12"));
//		genes.add(new Gene("8", "1", 0, 0, "p12"));
//		genes.add(new Gene("9", "1", 0, 0, "p12"));
//
//		HashSet<String> includedGenes = new HashSet<>();
//		for (int g = 1; g <= 8; ++g) {
//			includedGenes.add(String.valueOf(g));
//		}
//
//		HashMap<String, ArrayList<PathwayEnrichments.MetaGene>> metagenes = PathwayEnrichments.groupCorrelatedGenesPerChrArm(geneZscoresNullGwas, 0.9, genes, includedGenes, debugFolder, "NA", null);
//
//		Assert.assertTrue(metagenes.containsKey("1_q"));
//		Assert.assertTrue(metagenes.containsKey("1_p"));
//
//		Assert.assertEquals(metagenes.get("1_q").size(), 2);
//		Assert.assertEquals(metagenes.get("1_p").size(), 2);
//
//		Assert.assertTrue(metagenes.get("1_q").contains(new MetaGene(0,0,"1", "2", "4")));
//		Assert.assertTrue(metagenes.get("1_q").contains(new MetaGene(0,0,"3")));
//		Assert.assertTrue(metagenes.get("1_p").contains(new MetaGene(0,0,"5", "8")));
//		Assert.assertTrue(metagenes.get("1_p").contains(new MetaGene(0,0,"6", "7")));
//
//		//Meta genes are okay, now test collapse
//		DoubleMatrixDataset<String, String> mergedMean = PathwayEnrichments.collapseDatasetToMetaGenes(geneZscoresNullGwas, false, metagenes.values());
//
//		//view row to make sure order is the same
//		compareTwoMatrices(mergedMean, mergedMetaGenesMeanRef.viewRowSelection(mergedMean.getRowObjects()), 0.00001);
//		
//		
//		DoubleMatrixDataset<String, String> mergedZsum = PathwayEnrichments.collapseDatasetToMetaGenes(geneZscoresNullGwas, true, metagenes.values());
//
//		//view row to make sure order is the same
//		compareTwoMatrices(mergedZsum, mergedMetaGenesZscoreSumRef.viewRowSelection(mergedMean.getRowObjects()));
//		
//	}

	public static void compareTwoMatrices(DoubleMatrixDataset<String, String> m1, DoubleMatrixDataset<String, String> m2) {

		compareTwoMatrices(m1, m2, 0.00000001);

	}

	public static void compareTwoMatrices(DoubleMatrixDataset<String, String> m1, DoubleMatrixDataset<String, String> m2, double delta) {

		assertEquals(m1.rows(), m2.rows());
		assertEquals(m1.columns(), m2.columns());

		assertEquals(m1.getRowObjects(), m2.getRowObjects());
		assertEquals(m1.getColObjects(), m2.getColObjects());

		for (int r = 0; r < m1.rows(); ++r) {
			for (int c = 0; c < m1.columns(); ++c) {
				assertEquals(m1.getElementQuick(r, c), m2.getElementQuick(r, c), delta, "Difference at r: " + r + " c: " + c);
			}
		}

	}

}
