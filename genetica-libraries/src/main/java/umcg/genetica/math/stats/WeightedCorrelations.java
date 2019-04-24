/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class WeightedCorrelations {

	/**
	 * Weighted Pearson correlation
	 *
	 * Rows must be identical, only row count will be checked
	 *
	 * Columns in d1 will be columns in output, columns in d2 will be rows
	 *
	 * Based on code of Lude
	 *
	 * @param d1
	 * @param d2
	 * @param rowWeigths first column should contain weights
	 * @return
	 * @throws Exception
	 */
	public static DoubleMatrixDataset<String, String> weightedCorrelationColumnsOf2Datasets(DoubleMatrixDataset<String, String> d1, DoubleMatrixDataset<String, String> d2, DoubleMatrixDataset<String, String> rowWeigths) throws Exception {

		if (d1.rows() != d2.rows() || d1.rows() != rowWeigths.rows()) {
			throw new Exception("When correlating two datasets both should have identical number of rows. d1 has: " + d1.rows() + " d2 has: " + d2.rows() + " weights has: " + rowWeigths.rows() + " rows");
		}

		final DoubleMatrix2D d1Matrix = d1.getMatrix();
		final DoubleMatrix2D d2Matrix = d2.getMatrix();

		final DoubleMatrix1D weigths = rowWeigths.getCol(0);
		final double sumOfWeights = weigths.zSum();

		final DoubleMatrixDataset<String, String> correlations = new DoubleMatrixDataset<>(d2.getHashColsCopy(), d1.getHashColsCopy());

		final DoubleMatrix2D corMatrix = correlations.getMatrix();

		final int d1NrCols = d1.columns();
		final int d2NrCols = d2.columns();

		final int nrRows = d1.rows();

		final DoubleMatrix1D[] d1Cols = new DoubleMatrix1D[d1NrCols];
		final double[] d1ColsWeightedMean = new double[d1NrCols];
		for (int i = d1NrCols; --i >= 0;) {
			d1Cols[i] = d1Matrix.viewColumn(i);
			d1ColsWeightedMean[i] = weightedMean(d1Cols[i], weigths);
		}

		final DoubleMatrix1D[] d2Cols = new DoubleMatrix1D[d2NrCols];
		final double[] d2ColsWeightedMean = new double[d2NrCols];
		for (int i = d2NrCols; --i >= 0;) {
			d2Cols[i] = d2Matrix.viewColumn(i);
			d2ColsWeightedMean[i] = weightedMean(d2Cols[i], weigths);
		}

		for (int i = d1NrCols; --i >= 0;) {
			
			final DoubleMatrix1D col1 = d1Cols[i];
			final double wmX = d1ColsWeightedMean[i];
			
			for (int j = d2NrCols; --j >= 0;) {

				final DoubleMatrix1D col2 = d2Cols[j];

				final double wmY = d2ColsWeightedMean[j];

				double covXX = 0;
				double covXY = 0;
				double covYY = 0;

				for (int e = 0; e < nrRows; ++e) {
					covXX += weigths.getQuick(e) * (col1.getQuick(e) - wmX) * (col1.getQuick(e) - wmX);
					covXY += weigths.getQuick(e) * (col1.getQuick(e) - wmX) * (col2.getQuick(e) - wmY);
					covYY += weigths.getQuick(e) * (col2.getQuick(e) - wmY) * (col2.getQuick(e) - wmY);

				}

				covXX /= sumOfWeights;
				covXY /= sumOfWeights;
				covYY /= sumOfWeights;
				final double corr = covXY / (Math.sqrt(covXX * covYY));

				corMatrix.setQuick(j, i, corr);

			}
		}

		return correlations;

	}

	/**
	 * Weighted mean by Lude
	 *
	 * @param x
	 * @param weights
	 * @return
	 */
	public static double weightedMean(DoubleMatrix1D x, DoubleMatrix1D weights) {
		double m = 0;
		double sumWeights = 0;
		for (int s = 0; s < x.size(); s++) {
			m += x.getQuick(s) * weights.getQuick(s);
			sumWeights += weights.getQuick(s);
		}
		return m / sumWeights;
	}

}
