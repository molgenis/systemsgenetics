package deconvolution;

import java.io.IOException;

import org.apache.commons.lang3.NotImplementedException;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;

public class NotImplementedFunctions {
	/**
	 * Calculate the sum of squares, using Non-Negative Linear Regression, given a y expression vector with y ~
	 * model. This calculates the Beta parameters three times, for samples with GT = 0, GT = 1, GT = 2
	 * instead of using the interaction terms of the actual genotype.
	 * If no_intercept == true, remove the intercept (equivalent to y * ~ model -1 in R)
	 * 
	 * @param model An InteractionModel object including the y vector expression values and ObservedValues (model)
	 * Such that
	 * test_trait ~ geno_A + lymph% + geno_A:geno_B it can be for one QTL
	 * [[2, 43.4, 86.8], [2, 40.3, 80.6]], for another QTL [[0, 46.7, 0],
	 * [0, 51.5, 0] [0, 48.7, 0]] 
	 * 
	 * @return	The parameters per genotype dosage
	 */
	public static void calculateSumOfSquaresNNLS(LeastSquareModel model, Boolean plotBetaTimesVariables) throws IOException, IllegalAccessException {
		// OLS = Ordinary Least Squares
		OLSMultipleLinearRegression regression = new OLSMultipleLinearRegression();
		// if GetIntercept is false, remove the intercept (Beta1) from the linear model
		regression.setNoIntercept(model.getNoIntercept());
		regression.newSampleData(model.getExpessionValuesGt0(), model.getCellCountsGt0());
		regression.estimateRegressionParameters();
		throw new NotImplementedException("NNLS not implemneted");
	}
	public static double maximumLikelihoodEstimator() throws IOException, IllegalAccessException {
		/**
		 * Calculate the p-values of the betas by MLE estimation of parameters
		 */

		throw new NotImplementedException("mle not implemented yet");
	}
}
