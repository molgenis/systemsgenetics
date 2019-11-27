/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math;

import java.io.File;
import static org.testng.Assert.*;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class PcaColtNGTest {
	
	public PcaColtNGTest() {
	}

	@BeforeClass
	public static void setUpClass() throws Exception {
	}

	@BeforeMethod
	public void setUpMethod() throws Exception {
	}

	/**
	 * Test of getEigenvectors method, of class PcaColt.
	 */
	@Test
	public void testPca() throws Exception {
		File testMatrixFile = new File(this.getClass().getResource("/testMatrix.txt").toURI());
		
		DoubleMatrixDataset<String, String> testMatrix = DoubleMatrixDataset.loadDoubleTextData(testMatrixFile.getPath(), '\t');
		
		PcaColt pcaRes = new PcaColt(testMatrix, true);
		
		pcaRes.getEigenValues().printMatrix();
		
		pcaRes.getEigenvectors().printMatrix();
		
		pcaRes.getPcs().printMatrix();
		
	}
	
}
