package umcg.genetica.math.stats;

import org.apache.commons.math3.linear.SingularMatrixException;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import umcg.genetica.containers.Pair;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.util.*;
import java.util.concurrent.atomic.AtomicBoolean;

public class VIF {
	boolean debug = false;

	// assume covariates are on the columns!
	public DoubleMatrixDataset<String, String> vifCorrect(DoubleMatrixDataset<String, String> finalCovariates, double threshold) throws Exception {

		// no need to check for VIF when there is only one covariate
		if (finalCovariates.columns() <= 1) {
			return finalCovariates;
		}

//		System.out.println("VIF: " + finalCovariates.rows() + " x " + finalCovariates.columns());

		// determine variance inflation factor
//		System.out.println("Checking variance inflation factor...");
		HashSet<Integer> skipCol = new HashSet<>();
		boolean inflated = true;
		int iter = 0;

		// do some preqc
		for (int col = 0; col < finalCovariates.columns(); col++) {
			HashMap<Double, Integer> uniqueVals = new HashMap<>();
			double[] y = finalCovariates.getCol(col).toArray(); //[row];
			double var = Descriptives.variance(y);
			if (var == 0) {
				skipCol.add(col);
				if (debug) {
					System.out.println("Iteration: " + iter + "\tCovariate: " + finalCovariates.getColObjects().get(col) + "\thas zero variance.");
				}
			} else {
				// count unique values
				for (int d = 0; d < y.length; d++) {
					Integer q = uniqueVals.get(y[d]);
					if (q == null) {
						q = 0;
					}
					q++;
					uniqueVals.put(y[d], q);
				}

				if (uniqueVals.size() < 4) {
					// check whether each unique value is at least present 3 times
					AtomicBoolean b = new AtomicBoolean(true);

					DoubleMatrixDataset<String, String> finalCovariates1 = finalCovariates;
					int finalCol = col;
					int finalIter = iter;
					uniqueVals.forEach((k, v) -> {
						if (v < 3) {
							b.getAndSet(false);
							if (debug) {
								System.err.println("WARNING: Iteration: " + finalIter + "\tCovariate: " + finalCovariates1.getColObjects().get(finalCol) + "\thas " + uniqueVals.size() + " unique values, with " + k + " being present in " + v + " rows out of " + finalCovariates1.rows() + " ");
							}
						}
					});

					if (!b.get()) {
//						skipCol.add(col);
					}
				}
			}
		}


		iter = 1;
		finalCovariates = excludeCols(finalCovariates, skipCol);


		//		if (!skipCol.isEmpty()) {
//			System.out.println("There were problematic covariates. " + skipCol.size() + " covariates have been removed, " + (finalCovariates.columns()) + " remain.");
//		}

		// perform a pairwise correlation test first
		skipCol = new HashSet<>();
		for (int c1 = 0; c1 < finalCovariates.columns(); c1++) {
			if (!skipCol.contains(c1)) {
				double[] y = finalCovariates.getCol(c1).toArray(); //[row];
				for (int c2 = c1 + 1; c2 < finalCovariates.columns(); c2++) {
					double[] y2 = finalCovariates.getCol(c2).toArray(); //[row];
					double r = Correlation.correlate(y, y2);
					double rsq = r * r;
					if (rsq >= threshold) {
						if (debug) {
							System.out.println(finalCovariates.getColObjects().get(c1) + " correlates with " + finalCovariates.getColObjects().get(c2) + " - " + r + " removing latter");
						}
						skipCol.add(c2);
					}
				}
			}
		}

		finalCovariates = excludeCols(finalCovariates, skipCol);

		while (inflated) {
			skipCol = new HashSet<>();

//			ArrayList<Pair<Integer, Double>> colVarPairs = new ArrayList<Pair<Integer, Double>>();
//			for (int col = 0; col < finalCovariates.columns(); col++) {
//				HashMap<Double, Integer> uniqueVals = new HashMap<>();
//				double[] y = finalCovariates.getCol(col).toArray(); //[row];
//				double var = Descriptives.variance(y);
//				colVarPairs.add(new Pair<>(col, var, Pair.SORTBY.RIGHT));
//			}
//			Collections.sort(colVarPairs, Collections.reverseOrder());

//			for (Pair<Integer, Double> p : colVarPairs) {
//				System.out.println(p.getLeft() + "\t" + finalCovariates.getColObjects().get(p.getLeft()) + "\t" + p.getRight());
//			}

			for (int col = 0; col < finalCovariates.columns(); col++) {
//			for (int q = 0; q < colVarPairs.size(); q++) {
//				int col = colVarPairs.get(q).getLeft();
				try {
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

					double rsq = 1;

					rsq = ols.calculateRSquared();
					double vif = 1 / (1 - rsq);
					boolean alias = false;

					if (rsq > threshold) {
						alias = true;
						skipCol.add(col);
						if (debug) {
							System.out.println("VIF check: Iteration: " + iter + "\tCovariate: " + finalCovariates.getColObjects().get(col) + "\tRSq: " + rsq + "\tVIF: " + vif + "\tAliased: " + alias);
						}
						break;
					} else {
//						System.out.println("Iteration: " + iter + "\tCovariate: " + finalCovariates.getColObjects().get(col) + "\tRSq: " + rsq + "\tVIF: " + vif + "\tAliased: " + alias);
					}
				} catch (SingularMatrixException e) {
					if (debug) {
						System.out.println("Iteration: " + iter + "\tVIF correction produces singular matrix, when evaluating: " + finalCovariates.getColObjects().get(col) + "\tIgnoring covariate for now.");
					}
					// remove lowest variance covariate
					skipCol.add(col);
					break;
				}
			}

			if (skipCol.isEmpty()) {
				// System.out.println("There are no more collinear covariates. " + skipCol.size() + " covariates will be removed, " + (finalCovariates.columns() - skipCol.size()) + " remain.");
				inflated = false;
			} else {

				finalCovariates = excludeCols(finalCovariates, skipCol);
				if (debug) {
					System.out.println("There are collinear covariates. " + (finalCovariates.columns()) + " covariates remain.");
				}
				inflated = true;
			}
			iter++;


		}

		return finalCovariates;
	}

	private DoubleMatrixDataset<String, String> excludeCols
			(DoubleMatrixDataset<String, String> finalCovariates, HashSet<Integer> skipCol) throws Exception {
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

	private DoubleMatrixDataset<String, String> excludeRows
			(DoubleMatrixDataset<String, String> finalCovariates, HashSet<Integer> skipRow) throws Exception {
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

	public void setDebug(boolean b) {
		this.debug = b;
	}
}
