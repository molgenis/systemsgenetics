package umcg.genetica.math.stats;

import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.util.ArrayList;
import java.util.HashSet;

public class VIF {

	public DoubleMatrixDataset<String, String> vifCorrect(DoubleMatrixDataset<String, String> finalCovariates, double threshold) throws Exception {

		// determine variance inflation factor
		System.out.println("Checking variance inflation factor...");
		HashSet<Integer> skipRow = new HashSet<>();
		boolean inflated = true;
		int iter = 0;
		while (inflated) {
			skipRow = new HashSet<>();
			for (int row = 0; row < finalCovariates.rows(); row++) {
				OLSMultipleLinearRegression ols = new OLSMultipleLinearRegression();
				double[] y = finalCovariates.getRow(row).toArray(); //[row];

				// check if variance is >0
				double[][] otherCovariates = new double[finalCovariates.columns()][finalCovariates.rows() - 1];
				int rowctr = 0;

				for (int row2 = 0; row2 < finalCovariates.rows(); row2++) {
					if (row != row2) {
						for (int s = 0; s < finalCovariates.columns(); s++) {
							otherCovariates[s][rowctr] = finalCovariates.getElementQuick(row2, s);
						}
						rowctr++;
					}
				}

				ols.newSampleData(y, otherCovariates);

				double rsq = ols.calculateRSquared();
				double vif = 1 / (1 - rsq);
				boolean alias = false;

				if (rsq > threshold) {
					alias = true;
					skipRow.add(row);
					System.out.println("Iteration: " + iter + "\tCovariate: " + finalCovariates.getRowObjects().get(row) + "\tRSq: " + rsq + "\tVIF: " + vif + "\tAliased: " + alias);
					break;
				} else {
					System.out.println("Iteration: " + iter + "\tCovariate: " + finalCovariates.getRowObjects().get(row) + "\tRSq: " + rsq + "\tVIF: " + vif + "\tAliased: " + alias);
				}
			}

			if (skipRow.isEmpty()) {
				System.out.println("There are no more collinear covariates.");
				inflated = false;
			} else {
				finalCovariates = excludeRows(finalCovariates, skipRow);
				inflated = true;
			}
			iter++;
		}

		return finalCovariates;
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
