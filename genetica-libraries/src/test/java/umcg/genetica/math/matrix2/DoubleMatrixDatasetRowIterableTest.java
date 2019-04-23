/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.matrix2;

import java.io.File;
import java.net.URISyntaxException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Iterator;
import static org.testng.Assert.*;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;
import static umcg.genetica.math.matrix2.DoubleMatrixDatasetTest.compareTwoMatrices;

/**
 *
 * @author patri
 */
public class DoubleMatrixDatasetRowIterableTest {

	private File tmpOutputFolder;

	public DoubleMatrixDatasetRowIterableTest() {
	}

	@BeforeClass
	public static void setUpClass() throws Exception {
	}

	@BeforeMethod
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
	 * Test of getRows method, of class DoubleMatrixDatasetRowIterable.
	 */
	@Test
	public void test() throws URISyntaxException, Exception {

		File testMatrixFile = new File(this.getClass().getResource("/testMatrix.txt").toURI());
		
		DoubleMatrixDataset<String, String> testMatrix = DoubleMatrixDataset.loadDoubleTextData(testMatrixFile.getPath(), '\t');
		
		testMatrix.saveBinary(tmpOutputFolder.getAbsolutePath() + ".testMatrixBin");
		
		DoubleMatrixDatasetRowIterable matrixIterable = new DoubleMatrixDatasetRowIterable(tmpOutputFolder.getAbsolutePath() + ".testMatrixBin");
		
		assertEquals(matrixIterable.getNrRows(), testMatrix.rows());
		assertEquals(matrixIterable.getNrCols(), testMatrix.columns());
		
		double[][] dataLoaded = new double[matrixIterable.getNrRows()][0];
		
		Iterator<double[]> matrixIterator = matrixIterable.iterator();
		
		int row = 0;
		while(matrixIterator.hasNext()){
			
			dataLoaded[row++] = matrixIterator.next();
			
		}
	
		DoubleMatrixDataset<String, String> matrixFromIterator = new DoubleMatrixDataset<>(dataLoaded, new ArrayList(matrixIterable.getRows()), new ArrayList(matrixIterable.getCols()));
		
		compareTwoMatrices(matrixFromIterator, testMatrix);

	}

}
