/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.matrix2;

import java.io.IOException;
import java.util.LinkedHashMap;
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
public class DoubleMatrixDatasetRowCompressedReaderTest {
	
	public DoubleMatrixDatasetRowCompressedReaderTest() {
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

	@Test
	public void testSomeMethod() throws IOException, Exception {
		
		
		DoubleMatrixDataset<String, String> dummyData = new DoubleMatrixDataset<>(2, 3);
		dummyData.hashCols.put("C1", 0);
		dummyData.hashCols.put("C2", 1);
		dummyData.hashCols.put("C3", 2);
		dummyData.hashRows.put("R1", 0);
		dummyData.hashRows.put("R2", 1);
		dummyData.setElementQuick(0, 0, 1);
		dummyData.setElementQuick(0, 1, 2);
		dummyData.setElementQuick(0, 2, 3);
		dummyData.setElementQuick(1, 0, 4);
		dummyData.setElementQuick(1, 1, 5);
		dummyData.setElementQuick(1, 2, 6);
		
		dummyData.printMatrix();
		
		DoubleMatrixDatasetRowCompressedWriter.saveDataset("D:\\UMCG\\Genetica\\Projects\\tmp\\dummy", dummyData, "dummy", "rows", "cols");
		
			
		DoubleMatrixDataset<String, String> dummyData2 = new DoubleMatrixDatasetRowCompressedReader("D:\\UMCG\\Genetica\\Projects\\tmp\\dummy").loadFullDataset();
		dummyData.printMatrix();
		

		
		
		DoubleMatrixDatasetRowCompressedReader reader = new DoubleMatrixDatasetRowCompressedReader("D:\\UMCG\\Genetica\\Projects\\tmp\\test7.datg");

		long b = System.currentTimeMillis();
		DoubleMatrixDataset<String, String> z = reader.loadFullDataset();
		System.out.println(System.currentTimeMillis() - b);
			
		long a = System.currentTimeMillis();
		DoubleMatrixDataset<String, String> y = DoubleMatrixDataset.loadDoubleBinaryData("D:\\UMCG\\Genetica\\Projects\\tmp\\reactome_2020_07_18_raw");
		System.out.println(System.currentTimeMillis() - a);
		
	
		DoubleMatrixDatasetTest.compareTwoMatrices(y, z);
		
	}
	
}
