/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.matrix2;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashSet;
import static org.testng.Assert.*;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

/**
 *
 * @author patri
 */
public class DoubleMatrixDatasetTest {

	private File tmpOutputFolder;
	private File testMatrixFile;

	public DoubleMatrixDatasetTest() {
	}

	@BeforeClass
	public static void setUpClass() throws Exception {

	}

	@BeforeMethod
	public void setUpMethod() throws Exception {

		testMatrixFile = new File(this.getClass().getResource("/testMatrix.txt").toURI());

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
	 * Test of createRowForceNormalDuplicate method, of class
	 * DoubleMatrixDataset.
	 */
	@Test
	public void testCreateRowForceNormalDuplicate() {
		System.out.println("testCreateRowForceNormalDuplicate");
		ArrayList<String> rows = new ArrayList<>();
		ArrayList<String> cols = new ArrayList<>();

		rows.add("row1");
		rows.add("row2");
		rows.add("row3");
		rows.add("row4");

		cols.add("col1");
		cols.add("col2");
		cols.add("col3");
		cols.add("col4");
		cols.add("col5");

		DoubleMatrixDataset dataset = new DoubleMatrixDataset(rows, cols);

		dataset.setElementQuick(0, 0, 1);
		dataset.setElementQuick(0, 1, 2);
		dataset.setElementQuick(0, 2, 3);
		dataset.setElementQuick(0, 3, 4);
		dataset.setElementQuick(0, 4, 5);

		DoubleMatrixDataset datasetForNormal = dataset.createRowForceNormalDuplicate();

//		System.out.println(datasetForNormal.getElementQuick(0, 0));
//		System.out.println(datasetForNormal.getElementQuick(0, 1));
//		System.out.println(datasetForNormal.getElementQuick(0, 2));
//		System.out.println(datasetForNormal.getElementQuick(0, 3));
//		System.out.println(datasetForNormal.getElementQuick(0, 4));
//		
		assertEquals(datasetForNormal.getElementQuick(0, 0), 0.9736890569622489, 0.000001);
		assertEquals(datasetForNormal.getElementQuick(0, 4), 5.0263109430377515, 0.000001);

	}

	@Test
	public void testCalculateCorrelationMatrix() throws URISyntaxException, Exception {

		File testMatrixCorFile = new File(this.getClass().getResource("/testMatrixColumnCorMatrix.txt").toURI());

		DoubleMatrixDataset<String, String> testMatrix = DoubleMatrixDataset.loadDoubleTextData(testMatrixFile.getPath(), '\t');
		DoubleMatrixDataset<String, String> testMatrixRealCor = DoubleMatrixDataset.loadDoubleTextData(testMatrixCorFile.getPath(), '\t');

		DoubleMatrixDataset<String, String> corMatrix = testMatrix.calculateCorrelationMatrix();

		compareTwoMatrices(corMatrix, testMatrixRealCor);

	}

	@Test
	public void testCalculateCorrelationMatrixOnNormalizedColumns() throws URISyntaxException, Exception {

		File testMatrixCorFile = new File(this.getClass().getResource("/testMatrixColumnCorMatrix.txt").toURI());

		DoubleMatrixDataset<String, String> testMatrix = DoubleMatrixDataset.loadDoubleTextData(testMatrixFile.getPath(), '\t');

		testMatrix.normalizeColumns();

		DoubleMatrixDataset<String, String> testMatrixRealCor = DoubleMatrixDataset.loadDoubleTextData(testMatrixCorFile.getPath(), '\t');

		DoubleMatrixDataset<String, String> corMatrix = testMatrix.calculateCorrelationMatrixOnNormalizedColumns();

//		System.out.println("Calculated");
//		corMatrix.printMatrix();
//
//		System.out.println("");
//		System.out.println("Reference");
//		testMatrixRealCor.printMatrix();

		compareTwoMatrices(corMatrix, testMatrixRealCor);

	}

	@Test
	public void testCalculateCovarianceMatrix() throws Exception {

		File testMatrixCovFile = new File(this.getClass().getResource("/testMatrixColumnCovMatrix.txt").toURI());

		DoubleMatrixDataset<String, String> testMatrix = DoubleMatrixDataset.loadDoubleTextData(testMatrixFile.getPath(), '\t');
		System.out.println(testMatrixCovFile.getPath());
		DoubleMatrixDataset<String, String> testMatrixRealCov = DoubleMatrixDataset.loadDoubleTextData(testMatrixCovFile.getPath(), '\t');

		DoubleMatrixDataset<String, String> covMatrix = testMatrix.calculateCovarianceMatrix();
		//DoubleMatrixDataset<String, String> covMatrix2 = testMatrix.calculateCorrelationMatrixOnNormalizedColumns();

//		System.out.println("Calculated");
//		covMatrix.printMatrix();
//		System.out.println("");
//		System.out.println("Reference");
//		testMatrixRealCov.printMatrix();

//		System.out.println("Calculated2");
//		covMatrix2.printMatrix();
//		
		compareTwoMatrices(covMatrix, testMatrixRealCov, 0.5);

	}
	
	@Test
	public void correlateColumnsOf2ColumnNormalizedDatasets() throws Exception{
		
		DoubleMatrixDataset<String, String> testMatrix = DoubleMatrixDataset.loadDoubleTextData(testMatrixFile.getPath(), '\t');
		
		File testMatrixCorFile = new File(this.getClass().getResource("/testMatrixColumnCorMatrix.txt").toURI());
		
		testMatrix.normalizeColumns();

		DoubleMatrixDataset<String, String> testMatrixRealCor = DoubleMatrixDataset.loadDoubleTextData(testMatrixCorFile.getPath(), '\t');

		DoubleMatrixDataset<String, String> corMatrix = DoubleMatrixDataset.correlateColumnsOf2ColumnNormalizedDatasets(testMatrix, testMatrix);

//		System.out.println("Calculated");
//		corMatrix.printMatrix();
//
//		System.out.println("");
//		System.out.println("Reference");
//		testMatrixRealCor.printMatrix();

		compareTwoMatrices(corMatrix, testMatrixRealCor);

		
	}

	@Test
	public void testNormalizeRows() throws URISyntaxException, Exception {

		DoubleMatrixDataset<String, String> testMatrix = DoubleMatrixDataset.loadDoubleTextData(testMatrixFile.getPath(), '\t');

		DoubleMatrixDataset<String, String> testMatrixT = testMatrix.viewDice();

		testMatrixT.normalizeRows();

		File testMatrixScaleFile = new File(this.getClass().getResource("/testMatrixColumnScaledMatrix.txt").toURI());
		DoubleMatrixDataset<String, String> testMatrixRealScale = DoubleMatrixDataset.loadDoubleTextData(testMatrixScaleFile.getPath(), '\t');

		compareTwoMatrices(testMatrix, testMatrixRealScale);

	}

	@Test
	public void testNormalizeColumns() throws URISyntaxException, Exception {

		DoubleMatrixDataset<String, String> testMatrix = DoubleMatrixDataset.loadDoubleTextData(testMatrixFile.getPath(), '\t');

		testMatrix.normalizeColumns();

		File testMatrixScaleFile = new File(this.getClass().getResource("/testMatrixColumnScaledMatrix.txt").toURI());
		DoubleMatrixDataset<String, String> testMatrixRealScale = DoubleMatrixDataset.loadDoubleTextData(testMatrixScaleFile.getPath(), '\t');

		compareTwoMatrices(testMatrix, testMatrixRealScale);

	}

	@Test
	public void testSaveLoadBinaryMatrix() throws IOException, Exception {
		System.out.println("testSaveLoadBinaryMatrix");
		ArrayList<String> rows = new ArrayList<>();
		ArrayList<String> cols = new ArrayList<>();

		rows.add("row1");
		rows.add("row2");
		rows.add("row3");
		rows.add("row4");

		cols.add("col1");
		cols.add("col2");
		cols.add("col3");
		cols.add("col4");
		cols.add("col5");

		DoubleMatrixDataset dataset = new DoubleMatrixDataset(rows, cols);

		dataset.setElementQuick(0, 0, 1);
		dataset.setElementQuick(0, 1, 2);
		dataset.setElementQuick(0, 2, 3);
		dataset.setElementQuick(0, 3, 4);
		dataset.setElementQuick(0, 4, 5);

		dataset.setElementQuick(1, 3, 5.55);
		dataset.setElementQuick(2, 2, 6.66);
		dataset.setElementQuick(2, 3, -12.2);

		dataset.saveBinary(tmpOutputFolder.getAbsolutePath() + ".testBin");

		DoubleMatrixDataset<String, String> dataset2 = DoubleMatrixDataset.loadDoubleBinaryData(tmpOutputFolder.getAbsolutePath() + ".testBin");

		assertEquals(dataset2.rows(), 4);
		assertEquals(dataset2.columns(), 5);

		assertEquals(dataset2.getRowObjects().get(2), "row3");
		assertEquals(dataset2.getColObjects().get(4), "col5");

		assertEquals(dataset2.getElementQuick(0, 2), 3d);
		assertEquals(dataset2.getElementQuick(0, 1), 2d);
		assertEquals(dataset2.getElementQuick(2, 2), 6.66d);

		DoubleMatrixDataset<String, String> dataset3 = DoubleMatrixDataset.loadSubsetOfRowsBinaryDoubleData(tmpOutputFolder.getAbsolutePath() + ".testBin", new String[]{"row3", "row2"});

		assertEquals(dataset3.rows(), 2);
		assertEquals(dataset3.columns(), 5);

		assertEquals(dataset3.getRowObjects().get(0), "row3");
		assertEquals(dataset3.getRowObjects().get(1), "row2");
		assertEquals(dataset3.getColObjects().get(2), "col3");

		assertEquals(dataset3.getElementQuick(0, 2), 6.66d);
		assertEquals(dataset3.getElementQuick(0, 3), -12.2d);
		assertEquals(dataset3.getElementQuick(1, 3), 5.55d);
		assertEquals(dataset3.getElementQuick(1, 2), 0d);

		DoubleMatrixDataset<String, String> dataset4 = dataset2.viewRowSelection(new String[]{"row3", "row2"});

		assertEquals(dataset4.rows(), 2);
		assertEquals(dataset4.columns(), 5);

		assertEquals(dataset4.getRowObjects().get(0), "row3");
		assertEquals(dataset4.getRowObjects().get(1), "row2");
		assertEquals(dataset4.getColObjects().get(2), "col3");

		assertEquals(dataset4.getElementQuick(0, 2), 6.66d);
		assertEquals(dataset4.getElementQuick(0, 3), -12.2d);
		assertEquals(dataset4.getElementQuick(1, 3), 5.55d);
		assertEquals(dataset4.getElementQuick(1, 2), 0d);

		DoubleMatrixDatasetFastSubsetLoader subsetLoader = new DoubleMatrixDatasetFastSubsetLoader(tmpOutputFolder.getAbsolutePath() + ".testBin");
		DoubleMatrixDataset<String, String> dataset5 = subsetLoader.loadSubsetOfRowsBinaryDoubleData(new String[]{"row3", "row2"});

		assertEquals(dataset5.rows(), 2);
		assertEquals(dataset5.columns(), 5);

		assertEquals(dataset5.getRowObjects().get(0), "row3");
		assertEquals(dataset5.getRowObjects().get(1), "row2");
		assertEquals(dataset5.getColObjects().get(2), "col3");

		assertEquals(dataset5.getElementQuick(0, 2), 6.66d);
		assertEquals(dataset5.getElementQuick(0, 3), -12.2d);
		assertEquals(dataset5.getElementQuick(1, 3), 5.55d);
		assertEquals(dataset5.getElementQuick(1, 2), 0d);

	}

	@Test
	public void testSaveLoadTextMatrix() throws IOException, Exception {
		System.out.println("testSaveLoadTextMatrix");
		ArrayList<String> rows = new ArrayList<>();
		ArrayList<String> cols = new ArrayList<>();

		rows.add("row1");
		rows.add("row2");
		rows.add("row3");
		rows.add("row4");

		cols.add("col1");
		cols.add("col2");
		cols.add("col3");
		cols.add("col4");
		cols.add("col5");

		DoubleMatrixDataset dataset = new DoubleMatrixDataset(rows, cols);

		dataset.setElementQuick(0, 0, 1);
		dataset.setElementQuick(0, 1, 2);
		dataset.setElementQuick(0, 2, 3);
		dataset.setElementQuick(0, 3, 4);
		dataset.setElementQuick(0, 4, 5);

		dataset.setElementQuick(1, 3, 5.55);
		dataset.setElementQuick(2, 2, 6.66);
		dataset.setElementQuick(2, 3, -12.2);

		dataset.save(tmpOutputFolder.getAbsolutePath() + ".testText.txt");

		DoubleMatrixDataset<String, String> dataset2 = DoubleMatrixDataset.loadDoubleTextData(tmpOutputFolder.getAbsolutePath() + ".testText.txt", '\t');

		assertEquals(dataset2.rows(), 4);
		assertEquals(dataset2.columns(), 5);

		assertEquals(dataset2.getRowObjects().get(2), "row3");
		assertEquals(dataset2.getColObjects().get(4), "col5");

		assertEquals(dataset2.getElementQuick(0, 2), 3d);
		assertEquals(dataset2.getElementQuick(0, 1), 2d);
		assertEquals(dataset2.getElementQuick(2, 2), 6.66d);

		HashSet<String> rowsToLoad = new HashSet<>();
		rowsToLoad.add("row3");
		rowsToLoad.add("row2");

		DoubleMatrixDataset<String, String> dataset3 = DoubleMatrixDataset.loadSubsetOfTextDoubleData(tmpOutputFolder.getAbsolutePath() + ".testText.txt", '\t', rowsToLoad, null);

		assertEquals(dataset3.rows(), 2);
		assertEquals(dataset3.columns(), 5);

		assertEquals(dataset3.getRowObjects().get(0), "row2");
		assertEquals(dataset3.getRowObjects().get(1), "row3");
		assertEquals(dataset3.getColObjects().get(2), "col3");

		assertEquals(dataset3.getElementQuick(1, 2), 6.66d);
		assertEquals(dataset3.getElementQuick(1, 3), -12.2d);
		assertEquals(dataset3.getElementQuick(0, 3), 5.55d);
		assertEquals(dataset3.getElementQuick(0, 2), 0d);

		DoubleMatrixDataset<String, String> dataset4 = dataset3.viewRowSelection(new String[]{"row3", "row2"});

		assertEquals(dataset4.rows(), 2);
		assertEquals(dataset4.columns(), 5);

		assertEquals(dataset4.getRowObjects().get(0), "row3");
		assertEquals(dataset4.getRowObjects().get(1), "row2");
		assertEquals(dataset4.getColObjects().get(2), "col3");

		assertEquals(dataset4.getElementQuick(0, 2), 6.66d);
		assertEquals(dataset4.getElementQuick(0, 3), -12.2d);
		assertEquals(dataset4.getElementQuick(1, 3), 5.55d);
		assertEquals(dataset4.getElementQuick(1, 2), 0d);

		HashSet<String> colsToLoad = new HashSet<>();
		colsToLoad.add("col2");
		colsToLoad.add("col4");

		DoubleMatrixDataset<String, String> dataset5 = DoubleMatrixDataset.loadSubsetOfTextDoubleData(tmpOutputFolder.getAbsolutePath() + ".testText.txt", '\t', rowsToLoad, colsToLoad);

		assertEquals(dataset5.rows(), 2);
		assertEquals(dataset5.columns(), 2);

		assertEquals(dataset5.getRowObjects().get(0), "row2");
		assertEquals(dataset5.getRowObjects().get(1), "row3");
		assertEquals(dataset5.getColObjects().get(0), "col2");
		assertEquals(dataset5.getColObjects().get(1), "col4");

		assertEquals(dataset5.getElementQuick(0, 0), 0d);
		assertEquals(dataset5.getElementQuick(0, 1), 5.55d);
		assertEquals(dataset5.getElementQuick(1, 0), 0d);
		assertEquals(dataset5.getElementQuick(1, 1), -12.2d);

		DoubleMatrixDataset<String, String> dataset6 = dataset3.viewColSelection(new String[]{"col2", "col4"});

		assertEquals(dataset6.rows(), 2);
		assertEquals(dataset6.columns(), 2);

		assertEquals(dataset6.getRowObjects().get(0), "row2");
		assertEquals(dataset6.getRowObjects().get(1), "row3");
		assertEquals(dataset6.getColObjects().get(0), "col2");
		assertEquals(dataset6.getColObjects().get(1), "col4");

		assertEquals(dataset6.getElementQuick(0, 0), 0d);
		assertEquals(dataset6.getElementQuick(0, 1), 5.55d);
		assertEquals(dataset6.getElementQuick(1, 0), 0d);
		assertEquals(dataset6.getElementQuick(1, 1), -12.2d);

		DoubleMatrixDataset<String, String> dataset7 = dataset3.viewColSelection(new String[]{"col4", "col2"});

		assertEquals(dataset7.rows(), 2);
		assertEquals(dataset7.columns(), 2);

		assertEquals(dataset7.getRowObjects().get(0), "row2");
		assertEquals(dataset7.getRowObjects().get(1), "row3");
		assertEquals(dataset7.getColObjects().get(0), "col4");
		assertEquals(dataset7.getColObjects().get(1), "col2");

		assertEquals(dataset7.getElementQuick(0, 1), 0d);
		assertEquals(dataset7.getElementQuick(0, 0), 5.55d);
		assertEquals(dataset7.getElementQuick(1, 1), 0d);
		assertEquals(dataset7.getElementQuick(1, 0), -12.2d);

		DoubleMatrixDataset<String, String> dataset8 = DoubleMatrixDataset.loadSubsetOfTextDoubleData(tmpOutputFolder.getAbsolutePath() + ".testText.txt", '\t', null, colsToLoad);

		assertEquals(dataset8.rows(), 4);
		assertEquals(dataset8.columns(), 2);

		assertEquals(dataset8.getRowObjects().get(1), "row2");
		assertEquals(dataset8.getRowObjects().get(2), "row3");
		assertEquals(dataset8.getColObjects().get(0), "col2");
		assertEquals(dataset8.getColObjects().get(1), "col4");

		assertEquals(dataset8.getElementQuick(1, 0), 0d);
		assertEquals(dataset8.getElementQuick(1, 1), 5.55d);
		assertEquals(dataset8.getElementQuick(2, 0), 0d);
		assertEquals(dataset8.getElementQuick(2, 1), -12.2d);

	}

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
	
	public static void compareTwoMatricesIgnoreNaN(DoubleMatrixDataset<String, String> m1, DoubleMatrixDataset<String, String> m2, double delta) {

		assertEquals(m1.rows(), m2.rows());
		assertEquals(m1.columns(), m2.columns());

		for (int r = 0; r < m1.rows(); ++r) {
			for (int c = 0; c < m1.columns(); ++c) {
				
				double e1 = m1.getElementQuick(r, c);
				double e2 = m2.getElementQuick(r, c);
				
				if(Double.isNaN(e1) || Double.isNaN(e2)){
					continue;
				}
				
				assertEquals(e1, e2, delta, "Difference at r: " + r + " c: " + c);
			}
		}

	}

}
