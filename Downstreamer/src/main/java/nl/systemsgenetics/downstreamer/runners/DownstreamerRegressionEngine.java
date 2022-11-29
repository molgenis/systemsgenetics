package nl.systemsgenetics.downstreamer.runners;

import cern.colt.function.tdouble.DoubleDoubleFunction;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import cern.colt.matrix.tdouble.algo.decomposition.DenseDoubleEigenvalueDecomposition;
import cern.jet.math.tdouble.DoubleFunctions;
import nl.systemsgenetics.downstreamer.Downstreamer;
import nl.systemsgenetics.downstreamer.runners.options.OptionsModeRegress;
import nl.systemsgenetics.downstreamer.summarystatistic.LinearRegressionResult;
import org.apache.log4j.Logger;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.util.*;

public class DownstreamerRegressionEngine {

    private static final Logger LOGGER = Logger.getLogger(DownstreamerRegressionEngine.class);

    // Vector functions
    public static final DoubleDoubleFunction minus = (a, b) -> a - b;
    public static final DoubleDoubleFunction plus = (a, b) -> a + b;

    public static void run(OptionsModeRegress options) throws Exception {

        DoubleMatrixDataset<String, String> X = DoubleMatrixDataset.loadDoubleTextData(options.getExplanatoryVariables().getPath(), '\t');
        DoubleMatrixDataset<String, String> Y = DoubleMatrixDataset.loadDoubleTextData(options.getResponseVariable().getPath(), '\t');

        DoubleMatrixDataset<String, String> covariates = null;
        if (options.hasCovariates()) {
            covariates = DoubleMatrixDataset.loadDoubleTextData(options.getCovariates().getPath(), '\t');
        }

        if (options.hasSigma()) {
            DoubleMatrixDataset<String, String> Sigma = DoubleMatrixDataset.loadDoubleTextData(options.getSigma().getPath(), '\t');
            performDownstreamerRegression(X, Y, Sigma, options.getPercentageOfVariance(), options);
        } else {
            DoubleMatrixDataset<String, String> U = DoubleMatrixDataset.loadDoubleTextData(options.getEigenvectors().getPath(), '\t');
            DoubleMatrixDataset<String, String> L = DoubleMatrixDataset.loadDoubleTextData(options.getEigenvalues().getPath(), '\t');

            performDownstreamerRegression(X, Y, U, L, options.getPercentageOfVariance(), options);
        }

    }

    public static LinearRegressionResult performDownstreamerRegression(DoubleMatrixDataset<String, String> X,
                                                                DoubleMatrixDataset<String, String> Y,
                                                                DoubleMatrixDataset<String, String> Sigma,
                                                                double percentageOfVariance,
                                                                OptionsModeRegress optionsModeRegress) throws Exception {

        int nColumns = Sigma.columns();

        LOGGER.info("[" + Downstreamer.DATE_TIME_FORMAT.format(new Date()) + "] Running eigen decomposition on Sigma");
        //eigenvalues are ordered in oposite order to R! so from small to large
        DenseDoubleEigenvalueDecomposition eigen = new DenseDoubleEigenvalueDecomposition(Sigma.getMatrix());

        LOGGER.info("[" + Downstreamer.DATE_TIME_FORMAT.format(new Date()) + "] Done running eigen decomposition on Sigma");

        LinkedHashMap<String, Integer> eigenNames = new LinkedHashMap<>();
        LinkedHashMap<String, Integer> eigenNamesInverted = new LinkedHashMap<>();
        for (int i=0; i < nColumns; i++) {
            eigenNames.put("V" + ((nColumns-1) - i), i);
            eigenNamesInverted.put("V" + i, i);
        }

        if (eigen.getImagEigenvalues().zSum() > 0) {
            throw new Exception("Non real eigenvalues. Should not happen?");
        }

        // Re order and convert to DoubleMatrixDataset
        DoubleMatrixDataset<String, String> U = new DoubleMatrixDataset<>(Sigma.getHashRows(), eigenNames);
        U.setMatrix(eigen.getV());
        U.reorderCols(eigenNamesInverted);

        DoubleMatrixDataset<String, String> L = new DoubleMatrixDataset<>(nColumns, 1);
        L.setHashRows(eigenNames);
        L.setHashCols(new LinkedHashMap<String, Integer>(){{put("diag", 0);}});
        L.setMatrix(eigen.getRealEigenvalues().reshape((int)eigen.getRealEigenvalues().size(), 1));
        L.reorderRows(eigenNamesInverted);

        U.save(optionsModeRegress.getOutputBasePath() + "_U.txt");
        L.save(optionsModeRegress.getOutputBasePath() + "_L.txt");

        return performDownstreamerRegression(X, Y, U, L, percentageOfVariance, optionsModeRegress);
    }


    public static LinearRegressionResult performDownstreamerRegression(DoubleMatrixDataset<String, String> X,
                                                                DoubleMatrixDataset<String, String> Y,
                                                                DoubleMatrixDataset<String, String> U,
                                                                DoubleMatrixDataset<String, String> L,
                                                                double percentageOfVariance,
                                                                OptionsModeRegress optionsModeRegress) {

        DoubleMatrix1D varPerEigen = cumulativeSum(L.viewCol(0), true);

        List<String> colsToMaintain = new ArrayList<>();

        for (int i=0; i< varPerEigen.size(); i++) {
            if (percentageOfVariance > varPerEigen.get(i) ) {
                colsToMaintain.add(U.getColObjects().get(i));
            }
        }

        int numVectorsToMaintain = colsToMaintain.size();
        LOGGER.info("Maintaining " + numVectorsToMaintain + " eigenvectors explaining " + varPerEigen.get(numVectorsToMaintain-1) + " % of the variance.");

        // Subset U and L on the number of eigenvectors to maintain
        U = U.viewColSelection(colsToMaintain);
        L = L.viewRowSelection(colsToMaintain);

        // Pre-transpose U as this can be re-used as well
        DoubleMatrix2D UHatT = transpose(U.getMatrix());
        DoubleMatrix1D LHatInv = L.viewCol(0).copy();
        LHatInv.assign(DoubleFunctions.inv);

        // Pre-compute Y^ as this can be re-used
        DoubleMatrix2D YHat = mult(UHatT, Y.getMatrix());

        LOGGER.info("[" + Downstreamer.DATE_TIME_FORMAT.format(new Date()) + "] Starting regression for " + X.columns() + " pathways.");

        // TODO: parallelize. Altough not sure how this would work with the colt parralelization
        for (int curPathway=0; curPathway < X.columns(); curPathway++) {
            double[] curRes = downstreamerRegressionPrecomp(X.viewCol(curPathway).reshape(X.rows(), 1),
                    YHat,
                    UHatT,
                    diag(LHatInv),
                    true);

            LOGGER.debug("debug point");
        }

        LOGGER.info("[" + Downstreamer.DATE_TIME_FORMAT.format(new Date()) + "] Done with regression for " + X.columns() + " pathways.");

        // TODO: collare results into LinearRegressionResult
        return null;
    }

    /**
     * Generic implementation of the downstreamer regression model where:
     *      t(UHat) * Y
     *      t(UHat)
     *      solve(LHat)
     * have been pre-computed to save time
     * @param X               Predictors
     * @param YHat            t(UHat) * Y
     * @param UHatT           t(UHat)
     * @param LHatInv         1/diag(LHat)
     * @param fitIntercept    Should an intercept term be fitted
     * @return double[] with beta and standard erros in the form  [beta_1, beta_2 ... beta_x, se_1, se_2 ... se_x]
     */
    public static double[] downstreamerRegressionPrecomp(DoubleMatrix2D X, DoubleMatrix2D YHat, DoubleMatrix2D UHatT, DoubleMatrix2D LHatInv, boolean fitIntercept) {

        // X.hat <- U.hat.t %*% X
        double[] results = new double[X.columns()*2];
        DoubleMatrix2D XHat = mult(UHatT, X);

        if (fitIntercept) {
            XHat = DoubleFactory2D.dense.appendColumns(DoubleFactory2D.dense.make(XHat.rows(), 1, 1), XHat);
            results = new double[(X.columns()*2) + 2];
        }

        //-----------------------------------------------------------------------------------------
        // beta
        //-----------------------------------------------------------------------------------------
        //  a    <- crossprod(Y.hat, L.hat.inv) %*% X.hat
        //  b    <- solve(crossprod(X.hat, L.hat.inv) %*% X.hat)
        //  beta <- a %*% b
        DoubleMatrix2D a = mult(crossprod(YHat, LHatInv), XHat);
        DoubleMatrix2D b = inverse(mult(crossprod(XHat, LHatInv), XHat));
        DoubleMatrix2D beta = mult(a, b);

        //-----------------------------------------------------------------------------------------
        // Residuals
        //-----------------------------------------------------------------------------------------
        // Predicted Y
        DoubleMatrix1D predicted = tcrossprod(XHat, beta).viewColumn(0);
        // Original Y
        DoubleMatrix1D residuals = YHat.viewColumn(0).copy();
        // Subtract predicted from original
        residuals.assign(predicted, minus);

        // Residual sum of squares
        double rss = mult(crossprod(residuals, LHatInv), residuals.reshape((int) residuals.size(), 1)).get(0,0);

        // Divide rss by the degrees of freedom
        double rssW = rss / (XHat.rows() - XHat.columns());

        //-----------------------------------------------------------------------------------------
        // Standard error
        //-----------------------------------------------------------------------------------------
        // se     <- sqrt(diag(rss.w * b))
        DoubleMatrix1D se = diag(b.assign(DoubleFunctions.mult(rssW))).assign(DoubleFunctions.sqrt);

        // Put the results as an array in the form
        // [beta_1, beta_2 ... beta_x, se_1, se_2 ... se_x]
        for (int i=0; i< results.length / 2; i++) {
            results[i] = beta.get(0, i);
        }

        for (int i=0; i < se.size(); i++) {
            results[i+(results.length/2)] = se.get(i);
        }

        return results;
    }

    private static DoubleMatrix1D cumulativeSum(DoubleMatrix1D A, boolean asPercentage) {

        double totalSum = 1;
        if (asPercentage) {
            totalSum = A.zSum();
        }
        DoubleMatrix1D result = A.like();
        double runningSum =0;
        for (int i = 0; i < A.size(); i++) {
            runningSum = runningSum + A.get(i);
            result.set(i, runningSum / totalSum);
        }

        return result;
    }

    // Below are just calls to DenseDoubleAlgebra.DEFAULT to clean up the code above.
    /**
     * Shorthand for A * t(B)
     * @param A A matrix
     * @param B A matrix
     * @return
     */
    private static DoubleMatrix2D tcrossprod(DoubleMatrix2D A, DoubleMatrix2D B) {
        return DenseDoubleAlgebra.DEFAULT.mult(A, transpose(B));
    }

    /**
     * Shorthand for t(A) * B
     * @param A A matrix
     * @param B A matrix
     * @return
     */
    private static DoubleMatrix2D crossprod(DoubleMatrix2D A, DoubleMatrix2D B) {
        return DenseDoubleAlgebra.DEFAULT.mult(transpose(A), B);
    }

    /**
     * Shorthand for t(A) * B
     * @param A A column vector
     * @param B A matrix
     * @return
     */
    private static DoubleMatrix2D crossprod(DoubleMatrix1D A, DoubleMatrix2D B) {
        return DenseDoubleAlgebra.DEFAULT.mult(transpose( A.reshape((int) A.size(), 1)), B);
    }

    /**
     * Shorthand for A * B
     * @param A A matrix
     * @param B A matrix
     * @return
     */
    private static DoubleMatrix2D mult(DoubleMatrix2D A, DoubleMatrix2D B) {
        return DenseDoubleAlgebra.DEFAULT.mult(A, B);
    }

    /**
     * Shorthand for A^-1
     * @param A A matrix
     * @return
     */
    private static DoubleMatrix2D inverse(DoubleMatrix2D A) {
        return DenseDoubleAlgebra.DEFAULT.inverse(A);
    }


    /**
     * Shorthand for t(A)
     * @param A
     * @return
     */
    private static DoubleMatrix2D transpose(DoubleMatrix2D A) {
        return DenseDoubleAlgebra.DEFAULT.transpose(A);
    }

    /**
     * Shorthand for diag(A)
     * @param A
     * @return
     */
    private static DoubleMatrix1D diag(DoubleMatrix2D A) {
        return DoubleFactory2D.dense.diagonal(A);
    }

    /**
     * Shorthand for diag(A)
     * @param A
     * @return
     */
    private static DoubleMatrix2D diag(DoubleMatrix1D A) {
        return DoubleFactory2D.dense.diagonal(A);
    }


}
