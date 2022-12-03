package nl.systemsgenetics.downstreamer.runners;

import cern.colt.function.tdouble.DoubleDoubleFunction;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import cern.jet.math.tdouble.DoubleFunctions;
import me.tongfei.progressbar.ProgressBar;
import nl.systemsgenetics.downstreamer.Downstreamer;
import nl.systemsgenetics.downstreamer.io.IoUtils;
import nl.systemsgenetics.downstreamer.runners.options.OptionsModeRegress;
import nl.systemsgenetics.downstreamer.summarystatistic.LinearRegressionResult;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.log4j.Logger;


import org.ojalgo.function.constant.PrimitiveMath;
import org.ojalgo.matrix.Primitive64Matrix;
import org.ojalgo.matrix.QuaternionMatrix;
import org.ojalgo.matrix.store.SparseStore;
import org.ojalgo.structure.Mutate1D;
import org.ojalgo.structure.Mutate2D;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.util.*;

public class DownstreamerRegressionEngineOjAlgo {

    private static final Logger LOGGER = Logger.getLogger(DownstreamerRegressionEngineOjAlgo.class);

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

        DoubleMatrixDataset<String, String> U = DoubleMatrixDataset.loadDoubleData(options.getEigenvectors().getPath());
        DoubleMatrixDataset<String, String> L = DoubleMatrixDataset.loadDoubleData(options.getEigenvalues().getPath());

        if (rowsToSelect != null) {
            U = U.viewRowSelection(rowsToSelect);
        }

        LOGGER.info("Dim U:" + U.getMatrix().rows() + "x" + U.getMatrix().columns());
        LOGGER.info("Dim L:" + L.getMatrix().rows() + "x" + L.getMatrix().columns());

        LOGGER.info("Done loading datasets");


        output = performDownstreamerRegression(X, Y, covariates, U, L, options.getPercentageOfVariance(), false);


        // Save linear regression results
        output.save(options.getOutputBasePath(), false);
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


        LOGGER.info("Converting to ojAlgo");

        Primitive64Matrix nYhat = convertToPrimitive(YHat, false);
        Primitive64Matrix nUHatT = convertToPrimitive(UHatT, true);
        Primitive64Matrix nX = convertToPrimitive(X);
        Primitive64Matrix nLHatInv = diag(convertToPrimitive(LHatInv, false));

        LOGGER.info("Done converting to ojAlgo");
        ProgressBar pb = new ProgressBar("Linear regressions",  X.columns());
        for (int curPathway = 0; curPathway < X.columns(); curPathway++) {

            Primitive64Matrix XCur = nX.column(curPathway);
            double[] curRes = downstreamerRegressionPrecomp(XCur,
                    nYhat,
                    nUHatT,
                    nLHatInv,
                    fitIntercept);

            result.appendBetas(curPathway, ArrayUtils.subarray(curRes, 0, curRes.length/2));
            result.appendSe(curPathway, ArrayUtils.subarray(curRes, curRes.length/2, curRes.length));

            pb.step();
        }

        pb.close();
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
    public static double[] downstreamerRegressionPrecomp(Primitive64Matrix X, Primitive64Matrix YHat, Primitive64Matrix UHatT, Primitive64Matrix LHatInv, boolean fitIntercept) {

        // X.hat <- U.hat.t %*% X
        double[] results = new double[(int)X.countColumns() * 2];
        Primitive64Matrix XHat = mult(UHatT, X);

        //-----------------------------------------------------------------------------------------
        // beta
        //-----------------------------------------------------------------------------------------
        //  a    <- crossprod(Y.hat, L.hat.inv) %*% X.hat
        //  b    <- solve(crossprod(X.hat, L.hat.inv) %*% X.hat)
        //  beta <- a %*% b
        Primitive64Matrix a = mult(crossprod(YHat, LHatInv), XHat);
        Primitive64Matrix b = inverse(mult(crossprod(XHat, LHatInv), XHat));
        Primitive64Matrix beta = mult(a, b);

        //-----------------------------------------------------------------------------------------
        // Residuals
        //-----------------------------------------------------------------------------------------
        // Predicted Y
        Primitive64Matrix predicted = tcrossprod(XHat, beta);

        // Original Y
        Primitive64Matrix residuals = YHat.column(0);

        // Subtract predicted from original
        residuals = residuals.subtract(predicted);

        // Residual sum of squares
        double rss = mult(crossprod(residuals, LHatInv), residuals).get(0, 0);

        // Divide rss by the degrees of freedom
        double rssW = rss / (XHat.countRows() - XHat.countColumns());

        //-----------------------------------------------------------------------------------------
        // Standard error
        //-----------------------------------------------------------------------------------------
        // se     <- sqrt(diag(rss.w * b))
        Primitive64Matrix.DenseReceiver seTmp = diag(b.multiply(rssW)).copy();
        seTmp.modifyAll(PrimitiveMath.ROOT.parameter(2));
        Primitive64Matrix se = seTmp.get();
        //DoubleMatrix1D se = diag(b.assign(DoubleFunctions.mult(rssW))).assign(DoubleFunctions.sqrt);

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


    private static Primitive64Matrix convertToPrimitive(DoubleMatrix1D input, boolean sparse) {

        //PhysicalStore.Factory<Double, Primitive64Store> storeFactory = Primitive64Store.FACTORY;
        //Primitive64Store output = storeFactory.make(input.rows(), input.columns());

        if (sparse) {
           Primitive64Matrix.SparseReceiver output = Primitive64Matrix.FACTORY.makeSparse(input.size(), 1);

            for (int i=0; i < input.size(); i++) {
                output.set(i, input.get(i));
            }
            return output.get();
        } else  {Primitive64Matrix.DenseReceiver output = Primitive64Matrix.FACTORY.makeDense(input.size(), 1);
            for (int i=0; i < input.size(); i++) {
                output.set(i, input.get(i));
            }
            return output.get();
        }


    }


    private static Primitive64Matrix convertToPrimitive(DoubleMatrix2D input, boolean sparse) {

        //PhysicalStore.Factory<Double, Primitive64Store> storeFactory = Primitive64Store.FACTORY;
        //Primitive64Store output = storeFactory.make(input.rows(), input.columns());

        if (sparse) {
            Primitive64Matrix.SparseReceiver output = Primitive64Matrix.FACTORY.makeSparse(input.rows(), input.columns());

            for (int i=0; i < input.rows(); i++) {
                for (int j=0; j < input.columns(); j++) {
                    double val =  input.get(i, j);
                    if (val != 0) {
                        output.fillOne(i, j, val);
                    }
                    //double val = input.get(i, j);
                    //if (val != 0) {
                     //   output.set(i, j, val);
                    //}
                }
            }
            return output.get();
        } else  {
            Primitive64Matrix.DenseReceiver output = Primitive64Matrix.FACTORY.makeDense(input.rows(), input.columns());
            for (int i=0; i < input.rows(); i++) {
                for (int j=0; j < input.columns(); j++) {
                    output.set(i, j, input.get(i, j));
                }
            }

            return output.get();
        }

    }

    private static Primitive64Matrix convertToPrimitive(DoubleMatrixDataset<String, String> input) {

        //PhysicalStore.Factory<Double, Primitive64Store> storeFactory = Primitive64Store.FACTORY;
        //Primitive64Store output = storeFactory.make(input.rows(), input.columns());
        Primitive64Matrix.DenseReceiver output = Primitive64Matrix.FACTORY.makeDense(input.rows(), input.columns());

        for (int i=0; i < input.rows(); i++) {
            for (int j=0; j < input.columns(); j++) {
                output.set(i, j, input.getElementQuick(i, j));
            }
        }

        return output.get();

    }


    // Below are just calls to DenseDoubleAlgebra.DEFAULT to clean up the code above.

    /**
     * Shorthand for A * t(B)
     *
     * @param A A matrix
     * @param B A matrix
     * @return
     */
    private static Primitive64Matrix tcrossprod(Primitive64Matrix A, Primitive64Matrix B) {
        return A.multiply(B.transpose());
    }

    /**
     * Shorthand for t(A) * B
     *
     * @param A A matrix
     * @param B A matrix
     * @return
     */
    private static Primitive64Matrix crossprod(Primitive64Matrix A, Primitive64Matrix B) {
        return A.transpose().multiply(B);
    }

    /**
     * Shorthand for A * B
     *
     * @param A A matrix
     * @param B A matrix
     * @return
     */
    private static Primitive64Matrix mult(Primitive64Matrix A, Primitive64Matrix B) {
        return A.multiply(B);
    }

    /**
     * Shorthand for A^-1
     *
     * @param A A matrix
     * @return
     */
    private static Primitive64Matrix inverse(Primitive64Matrix A) {
        return A.invert();
    }


    /**
     * Shorthand for t(A)
     *
     * @param A
     * @return
     */
    private static Primitive64Matrix transpose(Primitive64Matrix A) {
        return A.transpose();
    }

    /**
     * Shorthand for diag(A)
     *
     * @param A
     * @return
     */
    private static Primitive64Matrix diag(Primitive64Matrix A) {
        Primitive64Matrix.SparseReceiver sparseReceiver = Primitive64Matrix.FACTORY.makeSparse(A.countRows(), A.countRows());
        sparseReceiver.fillDiagonal(A.columns(0));
        return sparseReceiver.get();
    }



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
