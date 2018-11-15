/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.matrix2;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import static org.testng.Assert.*;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

/**
 *
 * @author patri
 */
public class DoubleMatrixDatasetTest {

	public DoubleMatrixDatasetTest() {
	}

	@BeforeClass
	public static void setUpClass() throws Exception {
	}

	@BeforeMethod
	public void setUpMethod() throws Exception {
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

}
