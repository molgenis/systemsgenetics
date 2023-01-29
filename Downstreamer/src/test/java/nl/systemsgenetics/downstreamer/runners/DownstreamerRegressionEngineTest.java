/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.runners;

import cern.colt.function.tdouble.DoubleFunction;
import cern.jet.math.tdouble.DoubleFunctions;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import nl.systemsgenetics.downstreamer.gene.Gene;
import nl.systemsgenetics.downstreamer.io.BlockDiagonalDoubleMatrixProvider;
import nl.systemsgenetics.downstreamer.io.DoubleMatrixDatasetBlockDiagonalProvider;
import nl.systemsgenetics.downstreamer.io.IoUtils;
import static nl.systemsgenetics.downstreamer.runners.DownstreamerRegressionEngine.createBlockDiagonalIndexFromGenes2;
import nl.systemsgenetics.downstreamer.runners.options.OptionsModeRegress;
import static org.testng.Assert.*;
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
public class DownstreamerRegressionEngineTest {

	public DownstreamerRegressionEngineTest() {
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
	 * Test of blockDiagonalEigenDecomposition method, of class
	 * DownstreamerRegressionEngine.
	 */
	@Test
	public void testBlockDiagonalEigenDecomposition() throws Exception {

		File genesFile = new File(this.getClass().getResource("/random/genes.txt").toURI());
		final Map<String, List<Gene>> chrArmGeneMap = IoUtils.readGenesAsChrArmMap(genesFile);

		File correlationMatrixFile = new File(this.getClass().getResource("/random/random_true_correlation_matrix.txt.gz").toURI());
		DoubleMatrixDataset<String, String> correlationMatrix = DoubleMatrixDataset.loadDoubleData(correlationMatrixFile.getAbsolutePath());
		BlockDiagonalDoubleMatrixProvider provider = new DoubleMatrixDatasetBlockDiagonalProvider(correlationMatrix);

		ArrayList<String> allGenes = new ArrayList<>();
		for (List<Gene> armGenes : chrArmGeneMap.values()) {
			for (Gene gene : armGenes) {
				String geneName = gene.getGene();
				if (correlationMatrix.containsRow(geneName)) {
					allGenes.add(geneName);
				}
			}
		}
		//Do this suffle to make sure to test all the gene matchers and indexers
		Collections.shuffle(allGenes, new Random(42));

		LinkedHashMap<String, ArrayList<String>> index = createBlockDiagonalIndexFromGenes2(chrArmGeneMap, allGenes);

		File eigenvalues = new File(this.getClass().getResource("/random/genecor_eigenvalues.txt").toURI());
		DoubleMatrixDataset<String, String> expectedEigenvalues = DoubleMatrixDataset.loadDoubleData(eigenvalues.getAbsolutePath());

		File eigenvectors = new File(this.getClass().getResource("/random/genecor_eigenvectors.txt").toURI());
		//order using all genes arraylist
		DoubleMatrixDataset<String, String> expectedEigenVectors = DoubleMatrixDataset.loadDoubleData(eigenvectors.getAbsolutePath()).viewRowSelection(allGenes);

		boolean useJblas = true;

		DoubleMatrixDataset[] result = DownstreamerRegressionEngine.blockDiagonalEigenDecomposition(allGenes, provider, index, useJblas);

		DoubleMatrixDataset eigenVectors = result[1];

		columns:
		for (int c = 0; c < eigenVectors.columns(); ++c) {
			for (int r = 0; r < eigenVectors.rows(); ++r) {
				if (eigenVectors.getElementQuick(r, c) != 0) {
					if (eigenVectors.getElementQuick(r, c) * expectedEigenVectors.getElementQuick(r, c) < 0) {
						eigenVectors.getCol(c).assign(DoubleFunctions.neg);
					}
					continue columns;
				}
			}

		}

		compareTwoMatrices(result[0], expectedEigenvalues, 0.00000001);
		compareTwoMatrices(eigenVectors, expectedEigenVectors, 0.0000001);
		
		
		
//		useJblas = false;
//
//		result = DownstreamerRegressionEngine.blockDiagonalEigenDecomposition(allGenes, provider, index, useJblas);
//
//		eigenVectors = result[1];
//
//		columns:
//		for (int c = 0; c < eigenVectors.columns(); ++c) {
//			for (int r = 0; r < eigenVectors.rows(); ++r) {
//				if (eigenVectors.getElementQuick(r, c) != 0) {
//					if (eigenVectors.getElementQuick(r, c) * expectedEigenVectors.getElementQuick(r, c) < 0) {
//						eigenVectors.getCol(c).assign(DoubleFunctions.neg);
//					}
//					continue columns;
//				}
//			}
//
//		}
//
//		compareTwoMatrices(result[0], expectedEigenvalues, 0.1);
//		compareTwoMatrices(eigenVectors, expectedEigenVectors, 0.1);
//		

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
