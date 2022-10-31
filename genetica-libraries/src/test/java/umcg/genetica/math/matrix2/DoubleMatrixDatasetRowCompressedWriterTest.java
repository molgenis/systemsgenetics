/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.matrix2;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
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
public class DoubleMatrixDatasetRowCompressedWriterTest {
	
	public DoubleMatrixDatasetRowCompressedWriterTest() {
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
	 * Test of saveDataset method, of class DoubleMatrixDatasetRowCompressedWriter.
	 */
	@Test
	public void testSaveDataset() throws Exception {
		
		DoubleMatrixDataset<String, String> data = DoubleMatrixDataset.loadDoubleBinaryData("D:\\UMCG\\Genetica\\Projects\\tmp\\reactome_2020_07_18_raw");
		
		DoubleMatrixDatasetRowCompressedWriter.saveDataset("D:\\UMCG\\Genetica\\Projects\\tmp\\test4", data);
		
	}
	
}
