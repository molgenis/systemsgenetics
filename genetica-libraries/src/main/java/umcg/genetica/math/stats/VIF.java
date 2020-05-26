package umcg.genetica.math.stats;

import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.util.ArrayList;
import java.util.HashSet;

public class VIF {

	// assume covariates are on the columns!
	public DoubleMatrixDataset<String, String> vifCorrect(DoubleMatrixDataset<String, String> finalCovariates, double threshold) throws Exception {

		System.out.println("VIF: " + finalCovariates.rows() + " x " + finalCovariates.columns());

		// determine variance inflation factor
		System.out.println("Checking variance inflation factor...");
		HashSet<Integer> skipCol = new HashSet<>();
		boolean inflated = true;
		int iter = 0;
		while (inflated) {
			skipCol = new HashSet<>();

			for (int col = 0; col < finalCovariates.columns(); col++) {
				OLSMultipleLinearRegression ols = new OLSMultipleLinearRegression();
				double[] y = finalCovariates.getCol(col).toArray(); //[row];
				// check if variance is >0
				double[][] otherCovariates = new double[finalCovariates.rows()][finalCovariates.columns() - 1];
				int colctr = 0;
				for (int col2 = 0; col2 < finalCovariates.columns(); col2++) {
					if (col != col2) {
						for (int row = 0; row < finalCovariates.rows(); row++) {
							otherCovariates[row][colctr] = finalCovariates.getElementQuick(row, col2);
						}
						colctr++;
					}
				}
				ols.newSampleData(y, otherCovariates);

				double rsq = ols.calculateRSquared();
				double vif = 1 / (1 - rsq);
				boolean alias = false;

				if (rsq > threshold) {
					alias = true;
					skipCol.add(col);
					System.out.println("Iteration: " + iter + "\tCovariate: " + finalCovariates.getColObjects().get(col) + "\tRSq: " + rsq + "\tVIF: " + vif + "\tAliased: " + alias);
					break;
				} else {
					System.out.println("Iteration: " + iter + "\tCovariate: " + finalCovariates.getColObjects().get(col) + "\tRSq: " + rsq + "\tVIF: " + vif + "\tAliased: " + alias);
				}
			}

			if (skipCol.isEmpty()) {
				System.out.println("There are no more collinear covariates.");
				inflated = false;
			} else {
				finalCovariates = excludeCols(finalCovariates, skipCol);
				inflated = true;
			}
			iter++;


		}

		return finalCovariates;
	}

	private DoubleMatrixDataset<String, String> excludeCols(DoubleMatrixDataset<String, String> finalCovariates, HashSet<Integer> skipCol) throws Exception {
		DoubleMatrixDataset<String, String> tmp = new DoubleMatrixDataset<>(finalCovariates.rows(), finalCovariates.columns() - skipCol.size());
		tmp.setRowObjects(finalCovariates.getRowObjects());
		ArrayList<String> colObjs = new ArrayList<>();
		int colctr = 0;
		for (int col = 0; col < finalCovariates.columns(); col++) {
			if (!skipCol.contains(col)) {
				for (int row = 0; row < finalCovariates.rows(); row++) {
					tmp.setElementQuick(row, colctr, finalCovariates.getElementQuick(row, col));
				}
				colObjs.add(finalCovariates.getColObjects().get(col));
				colctr++;
			}
		}

		tmp.setColObjects(colObjs);
		return tmp;
	}

	private DoubleMatrixDataset<String, String> excludeRows(DoubleMatrixDataset<String, String> finalCovariates, HashSet<Integer> skipRow) throws Exception {
		// remove aliased row
		DoubleMatrixDataset<String, String> tmp = new DoubleMatrixDataset<>(finalCovariates.rows() - skipRow.size(), finalCovariates.columns());
		tmp.setColObjects(finalCovariates.getColObjects());
		ArrayList<String> rowObjs = new ArrayList<>();
		int rowctr = 0;
		for (int row = 0; row < finalCovariates.rows(); row++) {
			if (!skipRow.contains(row)) {
				for (int col = 0; col < finalCovariates.columns(); col++) {
					tmp.setElementQuick(rowctr, col, finalCovariates.getElementQuick(row, col));
				}
				rowObjs.add(finalCovariates.getRowObjects().get(row));
				rowctr++;
			}
		}

		tmp.setRowObjects(rowObjs);
		return tmp;
	}

}
