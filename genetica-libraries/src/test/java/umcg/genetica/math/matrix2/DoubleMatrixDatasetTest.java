/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.matrix2;

import java.io.File;
import java.io.IOException;
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

	public DoubleMatrixDatasetTest() {
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
	 * Test of createRowForceNormalDuplicate method, of class
	 * DoubleMatrixDataset.
	 */
	@Test
	public void testCreateRowForceNormalDuplicate() {

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
	public void testSaveLoadBinaryMatrix() throws IOException, Exception {

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

}
