/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats;

import java.io.File;
import static org.testng.Assert.*;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetTest;

/**
 *
 * @author patri
 */
public class PearsonRToZscoreBinnedNGTest {
	
	public PearsonRToZscoreBinnedNGTest() {
	}

	@BeforeClass
	public static void setUpClass() throws Exception {
	}

	@BeforeMethod
	public void setUpMethod() throws Exception {
	}


	/**
	 * Test of inplaceRToZ method, of class PearsonRToZscoreBinned.
	 * @throws java.lang.Exception
	 */
	@Test
	public void testInplaceRToZ() throws Exception {
		
		File testMatrixFile = new File(this.getClass().getResource("/testMatrix.txt").toURI());
		DoubleMatrixDataset<String, String> testMatrix = DoubleMatrixDataset.loadDoubleTextData(testMatrixFile.getPath(), '\t');
		
		File testMatrixRealCorZscoresFile = new File(this.getClass().getResource("/testMatrixColumnCorZscoreMatrix.txt").toURI());
		DoubleMatrixDataset<String, String> testMatrixRealCorZscores = DoubleMatrixDataset.loadDoubleTextData(testMatrixRealCorZscoresFile.getPath(), '\t');
		
		
		DoubleMatrixDataset<String, String> corMatrix = testMatrix.calculateCorrelationMatrix();
		
		PearsonRToZscoreBinned r2zScore = new PearsonRToZscoreBinned(10000, testMatrix.rows());
		
		r2zScore.inplaceRToZ(corMatrix);
		
		System.out.println("Real Z-score");
		testMatrixRealCorZscores.printMatrix();
		
		System.out.println("Calculated Z-score");
		corMatrix.printMatrix();
		
		DoubleMatrixDatasetTest.compareTwoMatricesIgnoreNaN(corMatrix, testMatrixRealCorZscores, 0.001);
		
	}
	
}
