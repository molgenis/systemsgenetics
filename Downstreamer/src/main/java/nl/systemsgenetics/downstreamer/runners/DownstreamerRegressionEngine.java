package nl.systemsgenetics.downstreamer.runners;

import cern.colt.function.tdouble.DoubleDoubleFunction;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import cern.colt.matrix.tdouble.algo.decomposition.DenseDoubleEigenvalueDecomposition;
import cern.jet.math.tdouble.DoubleFunctions;
import nl.systemsgenetics.downstreamer.Downstreamer;
import nl.systemsgenetics.downstreamer.io.IoUtils;
import nl.systemsgenetics.downstreamer.runners.options.OptionsModeRegress;
import nl.systemsgenetics.downstreamer.summarystatistic.LinearRegressionResult;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.log4j.Logger;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.util.*;

public class DownstreamerRegressionEngine {

    private static final Logger LOGGER = Logger.getLogger(DownstreamerRegressionEngine.class);

    // Vector functions
    public static final DoubleDoubleFunction minus = (a, b) -> a - b;
    public static final DoubleDoubleFunction plus = (a, b) -> a + b;

    public static void run(OptionsModeRegress options) throws Exception {

        LOGGER.info("Loading datasets");
        // TODO: currently doesnt just load the subset, so not the most efficient but easier to implemnt
        DoubleMatrixDataset<String, String> X = DoubleMatrixDataset.loadDoubleData(options.getExplanatoryVariables().getPath());
        DoubleMatrixDataset<String, String> Y = DoubleMatrixDataset.loadDoubleData(options.getResponseVariable().getPath());
        DoubleMatrixDataset<String, String> covariates = null;

        // Subset the data
        HashSet<String> rowsToSelect = null;
        HashSet<String> colsToSelect = null;

        if (options.hasRowIncludeFilter()) {
            LOGGER.info("Filtering rows to rows provided in -ro");
            rowsToSelect = new HashSet<>(IoUtils.readMatrixAnnotations(options.getRowIncludeFilter()));
            X = X.viewRowSelection(rowsToSelect);
            Y = Y.viewRowSelection(rowsToSelect);
        }

        if (options.hasColumnIncludeFilter()) {
            colsToSelect = new HashSet<>(IoUtils.readMatrixAnnotations(options.getColumnIncludeFilter()));
            X = X.viewColSelection(colsToSelect);
        }

        if (options.hasCovariates()) {
            covariates = DoubleMatrixDataset.loadDoubleData(options.getCovariates().getPath());
            if (rowsToSelect != null) {
                covariates = covariates.viewRowSelection(rowsToSelect);
            }
        }

        LOGGER.info("Dim X:" + X.getMatrix().rows() + "x" + X.getMatrix().columns());
        LOGGER.info("Dim Y:" + Y.getMatrix().rows() + "x" + Y.getMatrix().columns());

        // Depending on input, first decompse Sigma, or run with precomputed eigen decomp.
        LinearRegressionResult output;
        if (options.hasSigma()) {
            DoubleMatrixDataset<String, String> Sigma = DoubleMatrixDataset.loadDoubleData(options.getSigma().getPath());

            if (rowsToSelect != null) {
                Sigma = Sigma.viewSelection(rowsToSelect, rowsToSelect);
            }

            LOGGER.info("Dim Sigma:" + Sigma.getMatrix().rows() + "x" + Sigma.getMatrix().columns());

            LOGGER.info("Done loading datasets");
            output = performDownstreamerRegression(X, Y, covariates, Sigma, options.getPercentageOfVariance(),false, options);

        } else {
            DoubleMatrixDataset<String, String> U = DoubleMatrixDataset.loadDoubleData(options.getEigenvectors().getPath());
            DoubleMatrixDataset<String, String> L = DoubleMatrixDataset.loadDoubleData(options.getEigenvalues().getPath());

            if (rowsToSelect != null) {
                U = U.viewRowSelection(rowsToSelect);
            }

            LOGGER.info("Dim U:" + U.getMatrix().rows() + "x" + U.getMatrix().columns());
            LOGGER.info("Dim L:" + L.getMatrix().rows() + "x" + L.getMatrix().columns());

            LOGGER.info("Done loading datasets");
            output = performDownstreamerRegression(X, Y, covariates, U, L, options.getPercentageOfVariance(), true);
        }

        // Save linear regression results
        output.save(options.getOutputBasePath(), false);
    }

    /**
     * Run DS regression if only Sigma is given. I.e. run an eigen decomposition first.
     * @param X
     * @param Y
     * @param C
     * @param Sigma
     * @param percentageOfVariance
     * @param fitIntercept
     * @param optionsModeRegress
     * @return
     * @throws Exception
     */
    public static LinearRegressionResult performDownstreamerRegression(DoubleMatrixDataset<String, String> X,
                                                                       DoubleMatrixDataset<String, String> Y,
                                                                       DoubleMatrixDataset<String, String> C,
                                                                       DoubleMatrixDataset<String, String> Sigma,
                                                                       double percentageOfVariance,
                                                                       boolean fitIntercept,
                                                                       OptionsModeRegress optionsModeRegress) throws Exception {

        int nColumns = Sigma.columns();

        LOGGER.info("[" + Downstreamer.DATE_TIME_FORMAT.format(new Date()) + "] Running eigen decomposition on Sigma");
        //eigenvalues are ordered in oposite order to R! so from small to large
        DenseDoubleEigenvalueDecomposition eigen = new DenseDoubleEigenvalueDecomposition(Sigma.getMatrix());

        LOGGER.info("[" + Downstreamer.DATE_TIME_FORMAT.format(new Date()) + "] Done running eigen decomposition on Sigma");

        LinkedHashMap<String, Integer> eigenNames = new LinkedHashMap<>();
        LinkedHashMap<String, Integer> eigenNamesInverted = new LinkedHashMap<>();
        for (int i = 0; i < nColumns; i++) {
            eigenNames.put("V" + ((nColumns - 1) - i), i);
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
        L.setHashCols(new LinkedHashMap<String, Integer>() {{
            put("diag", 0);
        }});
        L.setMatrix(eigen.getRealEigenvalues().reshape((int) eigen.getRealEigenvalues().size(), 1));
        L.reorderRows(eigenNamesInverted);

        U.save(optionsModeRegress.getOutputBasePath() + "_U.txt");
        L.save(optionsModeRegress.getOutputBasePath() + "_L.txt");

        return performDownstreamerRegression(X, Y, C, U, L, percentageOfVariance, fitIntercept);
    }


    /**
     * Runs DS regression if eigendecompostion has been pre-computed. These are given by U and L. Expects that
     * these still need to be filtered. If they are pre-selected make sure to set percentageOfVariance to 1.
     * @param X
     * @param Y
     * @param C
     * @param U
     * @param L
     * @param percentageOfVariance
     * @param fitIntercept
     * @return
     */
    public static LinearRegressionResult performDownstreamerRegression(DoubleMatrixDataset<String, String> X,
                                                                       DoubleMatrixDataset<String, String> Y,
                                                                       DoubleMatrixDataset<String, String> C,
                                                                       DoubleMatrixDataset<String, String> U,
                                                                       DoubleMatrixDataset<String, String> L,
                                                                       double percentageOfVariance,
                                                                       boolean fitIntercept) {

        DoubleMatrix1D varPerEigen = cumulativeSum(L.viewCol(0), true);

        List<String> colsToMaintain = new ArrayList<>();

        for (int i = 0; i < varPerEigen.size(); i++) {
            if (percentageOfVariance > varPerEigen.get(i)) {
                colsToMaintain.add(U.getColObjects().get(i));
            }
        }

        int numVectorsToMaintain = colsToMaintain.size();
        LOGGER.info("Maintaining " + numVectorsToMaintain + " eigenvectors explaining " + varPerEigen.get(numVectorsToMaintain - 1) + " % of the variance.");

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

        // Determine the degrees of freedom. -1 for main_effect
        int degreesOfFreedom = numVectorsToMaintain - 1;

        // Set the names of the predictors in the correct order
        List<String> predictorNames = new ArrayList<>();
        if (fitIntercept) {
            predictorNames.add("intercept");
            degreesOfFreedom = degreesOfFreedom - 1;
        }
        predictorNames.add("main_effect");
        if (C != null) {
            predictorNames.addAll(C.getColObjects());
            degreesOfFreedom = degreesOfFreedom - C.columns();
        }

        LinearRegressionResult result = new LinearRegressionResult(X.getColObjects(), predictorNames, degreesOfFreedom);

        // TODO: parallelize. Altough not sure how this would work with the colt parralelization, might be redundant?
        for (int curPathway = 0; curPathway < X.columns(); curPathway++) {
            DoubleMatrix2D XCur = X.viewCol(curPathway).reshape(X.rows(), 1);

            // If covariates are provided add them
            if (C != null) {
                XCur = DoubleFactory2D.dense.appendColumns(XCur, C.getMatrix());
            }

            double[] curRes = downstreamerRegressionPrecomp(XCur,
                    YHat,
                    UHatT,
                    diag(LHatInv),
                    fitIntercept);

            result.appendBetas(curPathway, ArrayUtils.subarray(curRes, 0, curRes.length/2));
            result.appendSe(curPathway, ArrayUtils.subarray(curRes, curRes.length/2, curRes.length));

            LOGGER.debug("debug point");
        }

        LOGGER.info("[" + Downstreamer.DATE_TIME_FORMAT.format(new Date()) + "] Done with regression for " + X.columns() + " pathways.");

        return result;
    }

    /**
     * Generic implementation of the downstreamer regression model where:
     * t(UHat) * Y
     * t(UHat)
     * solve(LHat)
     * have been pre-computed to save time as these remain constant for all pathways in X.
     *
     * @param X            Predictors
     * @param YHat         t(UHat) * Y
     * @param UHatT        t(UHat)
     * @param LHatInv      1/diag(LHat)
     * @param fitIntercept Should an intercept term be fitted
     * @return double[] with beta and standard erros in the form  [beta_1, beta_2 ... beta_x, se_1, se_2 ... se_x]
     */
    public static double[] downstreamerRegressionPrecomp(DoubleMatrix2D X, DoubleMatrix2D YHat, DoubleMatrix2D UHatT, DoubleMatrix2D LHatInv, boolean fitIntercept) {

        // X.hat <- U.hat.t %*% X
        double[] results = new double[X.columns() * 2];
        DoubleMatrix2D XHat = mult(UHatT, X);

        if (fitIntercept) {
            XHat = DoubleFactory2D.dense.appendColumns(DoubleFactory2D.dense.make(XHat.rows(), 1, 1), XHat);
            results = new double[(X.columns() * 2) + 2];
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
        double rss = mult(crossprod(residuals, LHatInv), residuals.reshape((int) residuals.size(), 1)).get(0, 0);

        // Divide rss by the degrees of freedom
        double rssW = rss / (XHat.rows() - XHat.columns());

        //-----------------------------------------------------------------------------------------
        // Standard error
        //-----------------------------------------------------------------------------------------
        // se     <- sqrt(diag(rss.w * b))
        DoubleMatrix1D se = diag(b.assign(DoubleFunctions.mult(rssW))).assign(DoubleFunctions.sqrt);

        // Put the results as an array in the form
        // [beta_1, beta_2 ... beta_x, se_1, se_2 ... se_x]
        for (int i = 0; i < results.length / 2; i++) {
            results[i] = beta.get(0, i);
        }

        for (int i = 0; i < se.size(); i++) {
            results[i + (results.length / 2)] = se.get(i);
        }

        return results;
    }

    private static DoubleMatrix1D cumulativeSum(DoubleMatrix1D A, boolean asPercentage) {

        double totalSum = 1;
        if (asPercentage) {
            totalSum = A.zSum();
        }
        DoubleMatrix1D result = A.like();
        double runningSum = 0;
        for (int i = 0; i < A.size(); i++) {
            runningSum = runningSum + A.get(i);
            result.set(i, runningSum / totalSum);
        }

        return result;
    }

    // Below are just calls to DenseDoubleAlgebra.DEFAULT to clean up the code above.

    /**
     * Shorthand for A * t(B)
     *
     * @param A A matrix
     * @param B A matrix
     * @return
     */
    private static DoubleMatrix2D tcrossprod(DoubleMatrix2D A, DoubleMatrix2D B) {
        return DenseDoubleAlgebra.DEFAULT.mult(A, transpose(B));
    }

    /**
     * Shorthand for t(A) * B
     *
     * @param A A matrix
     * @param B A matrix
     * @return
     */
    private static DoubleMatrix2D crossprod(DoubleMatrix2D A, DoubleMatrix2D B) {
        return DenseDoubleAlgebra.DEFAULT.mult(transpose(A), B);
    }

    /**
     * Shorthand for t(A) * B
     *
     * @param A A column vector
     * @param B A matrix
     * @return
     */
    private static DoubleMatrix2D crossprod(DoubleMatrix1D A, DoubleMatrix2D B) {
        return DenseDoubleAlgebra.DEFAULT.mult(transpose(A.reshape((int) A.size(), 1)), B);
    }

    /**
     * Shorthand for A * B
     *
     * @param A A matrix
     * @param B A matrix
     * @return
     */
    private static DoubleMatrix2D mult(DoubleMatrix2D A, DoubleMatrix2D B) {
        return DenseDoubleAlgebra.DEFAULT.mult(A, B);
    }

    /**
     * Shorthand for A^-1
     *
     * @param A A matrix
     * @return
     */
    private static DoubleMatrix2D inverse(DoubleMatrix2D A) {
        return DenseDoubleAlgebra.DEFAULT.inverse(A);
    }


    /**
     * Shorthand for t(A)
     *
     * @param A
     * @return
     */
    private static DoubleMatrix2D transpose(DoubleMatrix2D A) {
        return DenseDoubleAlgebra.DEFAULT.transpose(A);
    }

    /**
     * Shorthand for diag(A)
     *
     * @param A
     * @return
     */
    private static DoubleMatrix1D diag(DoubleMatrix2D A) {
        return DoubleFactory2D.dense.diagonal(A);
    }

    /**
     * Shorthand for diag(A)
     *
     * @param A
     * @return
     */
    private static DoubleMatrix2D diag(DoubleMatrix1D A) {
        return DoubleFactory2D.dense.diagonal(A);
    }


}
