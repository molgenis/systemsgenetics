/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats;

import java.io.File;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import static umcg.genetica.math.matrix2.DoubleMatrixDatasetTest.compareTwoMatrices;

/**
 *
 * @author patri
 */
public class WeightedCorrelationsTest {

	public WeightedCorrelationsTest() {
	}

	@BeforeClass
	public static void setUpClass() throws Exception {
	}

	@BeforeMethod
	public void setUpMethod() throws Exception {
	}

	/**
	 * Test of weightedCorrelationColumnsOf2Datasets method, of class
	 * WeightedCorrelations.
	 */
	@Test
	public void testWeightedCorrelationColumnsOf2Datasets() throws Exception {

		File testMatrixFile = new File(this.getClass().getResource("/testMatrix.txt").toURI());

		DoubleMatrixDataset<String, String> testMatrix = DoubleMatrixDataset.loadDoubleTextData(testMatrixFile.getPath(), '\t');

		File testMatrixWeightedCorFile = new File(this.getClass().getResource("/testMatrixColumnWeightedCorMatrix.txt").toURI());
		DoubleMatrixDataset<String, String> testMatrixRealWCor = DoubleMatrixDataset.loadDoubleTextData(testMatrixWeightedCorFile.getPath(), '\t');

		File weightsFile = new File(this.getClass().getResource("/randomWeights.txt").toURI());
		DoubleMatrixDataset<String, String> weights = DoubleMatrixDataset.loadDoubleTextData(weightsFile.getPath(), '\t');

		DoubleMatrixDataset<String, String> testMatrixWCor = WeightedCorrelations.weightedCorrelationColumnsOf2Datasets(testMatrix, testMatrix, weights);

		System.out.println("Calculated");
		testMatrixWCor.printMatrix();

		System.out.println("");
		System.out.println("Reference");
		testMatrixRealWCor.printMatrix();
		
		
		compareTwoMatrices(testMatrixWCor, testMatrixRealWCor);

	}

}
