package nl.systemsgenetics.downstreamer.runners;

import cern.colt.function.tdouble.DoubleDoubleFunction;
import cern.colt.function.tdouble.DoubleProcedure;
import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import cern.jet.math.tdouble.DoubleFunctions;
import me.tongfei.progressbar.ProgressBar;
import nl.systemsgenetics.downstreamer.Downstreamer;
import nl.systemsgenetics.downstreamer.gene.Gene;
import nl.systemsgenetics.downstreamer.io.IoUtils;
import nl.systemsgenetics.downstreamer.runners.options.OptionsModeRegress;
import nl.systemsgenetics.downstreamer.summarystatistic.LinearRegressionResult;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.log4j.Logger;
import org.jblas.DoubleMatrix;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.util.*;

public class DownstreamerRegressionEnginePColt {

    private static final Logger LOGGER = Logger.getLogger(DownstreamerRegressionEnginePColt.class);

    // Vector functions
    public static final DoubleDoubleFunction minus = (a, b) -> a - b;
    public static final DoubleDoubleFunction plus = (a, b) -> a + b;

    public static void run(OptionsModeRegress options) throws Exception {

        LOGGER.warn("Assumes input data are in the same order unless -ro is specified. Order in -ro should match with -u and -l");
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
            throw new Exception("Covariates not yet supported");

/*            covariates = DoubleMatrixDataset.loadDoubleData(options.getCovariates().getPath());
            if (rowsToSelect != null) {
                covariates = covariates.viewRowSelection(rowsToSelect);
            }*/
        }

        List<int[]> blockDiagonalIndices = null;
        List<String> rowOrder = X.getRowObjects();

        if (options.hasGenes()) {

            // Determine the indices in the matrix that lie together on a chromosome arm
            blockDiagonalIndices = new ArrayList<>();
            Map<String, List<Integer>> blockDiagonalIndicesTmp = new HashMap<>();
            Map<String, List<Gene>> genes = IoUtils.readGenesAsChrArmMap(options.getGenes());

            for (String key : genes.keySet()) {
                for (Gene curGene : genes.get(key)) {
                    int idx = rowOrder.indexOf(curGene.getGene());
                    if (idx > -1) {
                        if (!blockDiagonalIndicesTmp.containsKey(key)) {
                            blockDiagonalIndicesTmp.put(key, new ArrayList<>());
                        }
                        blockDiagonalIndicesTmp.get(key).add(idx);
                    }
                }
            }

            // Convert to int[] so it's easier to use with DoubleMatrix2D;
            for (List<Integer> curList : blockDiagonalIndicesTmp.values()) {
                blockDiagonalIndices.add(curList.stream().mapToInt(i -> i).toArray());
            }

            LOGGER.info("Dim Y:" + Y.getMatrix().rows() + "x" + Y.getMatrix().columns());
        }

        LOGGER.info("Dim X:" + X.getMatrix().rows() + "x" + X.getMatrix().columns());
        LOGGER.info("Dim Y:" + Y.getMatrix().rows() + "x" + Y.getMatrix().columns());

        // Depending on input, first decompse Sigma, or run with precomputed eigen decomp.
        LinearRegressionResult output;
        if (options.hasSigma()) {
            throw new Exception("Not supported atm");
        } else {
            DoubleMatrixDataset<String, String> U = DoubleMatrixDataset.loadDoubleData(options.getEigenvectors().getPath());
            DoubleMatrixDataset<String, String> L = DoubleMatrixDataset.loadDoubleData(options.getEigenvalues().getPath());

            if (rowsToSelect != null) {
                U = U.viewRowSelection(rowsToSelect);
            }

            LOGGER.info("Dim U:" + U.getMatrix().rows() + "x" + U.getMatrix().columns());
            LOGGER.info("Dim L:" + L.getMatrix().rows() + "x" + L.getMatrix().columns());

            LOGGER.info("Done loading datasets");
            output = performDownstreamerRegression(X, Y, covariates, U, L, blockDiagonalIndices, options.getPercentageOfVariance(), true);
        }

        // Save linear regression results
        output.save(options.getOutputBasePath(), false);
    }


    /**
     * Runs DS regression if eigendecompostion has been pre-computed. These are given by U and L. Expects that
     * these still need to be filtered. If they are pre-selected make sure to set percentageOfVariance to 1.
     *
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
                                                                       List<int[]> blockDiagonalIndices,
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
        //UHatT = DoubleFactory2D.sparse.make(UHatT.toArray());

        DoubleMatrix1D LHatInv = L.viewCol(0).copy();
        LHatInv.assign(DoubleFunctions.inv);
        DoubleMatrix2D LHatInvM = diag(LHatInv);

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

        DoubleMatrix2D XhatCache = null;
        if (blockDiagonalIndices != null) {
            LOGGER.info("[" + Downstreamer.DATE_TIME_FORMAT.format(new Date()) + "] Precomputing XHat ");
            XhatCache = multBlockDiagonalSubsetPerColJblas(UHatT, X.getMatrix(), blockDiagonalIndices);
            LOGGER.info("[" + Downstreamer.DATE_TIME_FORMAT.format(new Date()) + "] Done ");
        }

        ProgressBar pb = new ProgressBar("Linear regressions", X.columns());
        for (int curPathway = 0; curPathway < X.columns(); curPathway++) {
            //DoubleMatrix2D XCur = X.viewCol(curPathway).reshape(X.rows(), 1);
            // DoubleMatrix2D XCur = X.viewColAsMmatrix(curPathway);
            DoubleMatrix2D XCur;

            if (XhatCache == null) {
                 XCur = X.viewColAsMmatrix(curPathway);

                // If covariates are provided add them
                if (C != null) {
                    XCur = DoubleFactory2D.dense.appendColumns(XCur, C.getMatrix());
                }

                XCur = mult(UHatT, XCur);
            } else {
                XCur = XhatCache.viewColumn(curPathway).reshape(XhatCache.rows(), 1);
            }

            double[] curRes = downstreamerRegressionPrecomp(XCur,
                    YHat,
                    UHatT,
                    LHatInvM,
                    blockDiagonalIndices,
                    fitIntercept);

            result.appendBetas(curPathway, ArrayUtils.subarray(curRes, 0, curRes.length / 2));
            result.appendSe(curPathway, ArrayUtils.subarray(curRes, curRes.length / 2, curRes.length));

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
     * @param XHat         Predictors
     * @param YHat         t(UHat) * Y
     * @param UHatT        t(UHat)
     * @param LHatInv      1/diag(LHat)
     * @param fitIntercept Should an intercept term be fitted
     * @return double[] with beta and standard erros in the form  [beta_1, beta_2 ... beta_x, se_1, se_2 ... se_x]
     */
    public static double[] downstreamerRegressionPrecomp(DoubleMatrix2D XHat, DoubleMatrix2D YHat, DoubleMatrix2D UHatT, DoubleMatrix2D LHatInv, Collection<int[]> blockDiagonalIndices, boolean fitIntercept) {

        double[] results = new double[XHat.columns() * 2];

        // if (blockDiagonalIndices != null) {
        //    XHat = multBlockDiagonalSubsetAOnCols(UHatT, XHat, blockDiagonalIndices);
        //} else {
        //    XHat = mult(UHatT, XHat);
        // }

        if (fitIntercept) {
            results = new double[(XHat.columns() * 2) + 2];
            XHat = DoubleFactory2D.dense.appendColumns(DoubleFactory2D.dense.make(XHat.rows(), 1, 1), XHat);
        }

        //-----------------------------------------------------------------------------------------
        // beta
        //-----------------------------------------------------------------------------------------
        //  a    <- crossprod(Y.hat, L.hat.inv) %*% X.hat
        //  b    <- solve(crossprod(X.hat, L.hat.inv) %*% X.hat)
        //  beta <- a %*% b

        //DoubleMatrix2D a = mult(crossprod(YHat, LHatInv), XHat);
        //DoubleMatrix2D b = inverse(mult(crossprod(XHat, LHatInv), XHat));

        // Alternative for YHat %*% LHatInv, bit quicker
        DoubleMatrix2D a = mult(transpose(prodDiag(YHat, diag(LHatInv))), XHat);
        DoubleMatrix2D b = inverse(mult(transpose(prodDiag(XHat, diag(LHatInv))), XHat));
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
        //double rss = mult(crossprod(residuals, LHatInv), residuals.reshape((int) residuals.size(), 1)).get(0, 0);

        // Alternative for residuals %*% LHatInv, bit quicker,
        double rss = mult(prodDiagT(residuals, diag(LHatInv)), residuals);
        //DoubleMatrix2D tmp = residuals.reshape((int) residuals.size(), 1);
        //double rss = mult(transpose(prodDiag(tmp, diag(LHatInv))), tmp).get(0,0);

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
        return mult(A, transpose(B));
    }

    /**
     * Shorthand for t(A) * B
     *
     * @param A A matrix
     * @param B A matrix
     * @return
     */
    private static DoubleMatrix2D crossprod(DoubleMatrix2D A, DoubleMatrix2D B) {
        return mult(transpose(A), B);
    }

    /**
     * Shorthand for t(A) * B
     *
     * @param A A column vector
     * @param B A matrix
     * @return
     */
    private static DoubleMatrix2D crossprod(DoubleMatrix1D A, DoubleMatrix2D B) {
        return mult(transpose(A.reshape((int) A.size(), 1)), B);
    }


    private static DoubleMatrix2D multBlockDiagonalSubsetPerColJblas(DoubleMatrix2D Ac, DoubleMatrix2D Bc, List<int[]> index) {

        DoubleMatrix2D C = DoubleFactory2D.dense.make(Ac.rows(), Bc.columns());

        DoubleMatrix A = new DoubleMatrix(Ac.toArray());
        DoubleMatrix B = new DoubleMatrix(Bc.toArray());

        //ProgressBar pb = new ProgressBar("Precomp XHat", index.size() * B.columns);
        ProgressBar pb = new ProgressBar("Precomp XHat", B.columns * index.size());

        DoubleMatrix[] Acache = new DoubleMatrix[index.size()];
        DoubleMatrix[] Bcache = new DoubleMatrix[index.size()];

        int[][] AcacheIndices = new int[index.size()][];

        for (int i = 0; i < index.size(); i++) {
            Acache[i] = A.getColumns(index.get(i));
            DoubleMatrix tmp = Acache[i].rowSums();
            ArrayList<Integer> tmpIndices = new ArrayList<>();
            for (int j = 0; j < tmp.rows; j++) {
                if (tmp.get(j, 0) != 0) {
                    tmpIndices.add(j);
                }
            }

            AcacheIndices[i] = tmpIndices.stream().mapToInt(p -> p).toArray();
            Acache[i] = Acache[i].getRows(AcacheIndices[i]);
            Bcache[i] = B.getRows(index.get(i));
        }

        for (int i = 0; i < index.size(); i++) {
            DoubleMatrix Ctmp = new DoubleMatrix(Acache[i].rows, 1);

            for (int col = 0; col < B.columns; col++) {

                Acache[i].mmuli(Bcache[i].getColumn(col), Ctmp);
                for (int j = 0; j < Ctmp.rows; j++) {
                    int idx = AcacheIndices[i][j];
                    C.setQuick(idx, col, C.getQuick(idx, col) + Ctmp.get(j, 0));
                }
                pb.step();
            }

        }
        pb.close();

        return C;
    }

    private static DoubleMatrix2D multBlockDiagonalSubsetPerCol(DoubleMatrix2D A, DoubleMatrix2D B, List<int[]> index) {

        DoubleMatrix2D C = DoubleFactory2D.dense.make(A.rows(), B.columns());
        DoubleMatrix2D tmp = DoubleFactory2D.dense.make(A.rows(), 1);

        //ProgressBar pb = new ProgressBar("Precomp XHat", index.size() * B.columns);
        ProgressBar pb = new ProgressBar("Precomp XHat", B.columns() * index.size());

        DoubleMatrix2D[] Acache = new DoubleMatrix2D[index.size()];
        DoubleMatrix2D[] Bcache = new DoubleMatrix2D[index.size()];

        int[][] AcacheIndices = new int[index.size()][];

        for (int i = 0; i < index.size(); i++) {
            Acache[i] = A.viewSelection(null, index.get(i));
            ArrayList<Integer> tmpIndices = new ArrayList<>();

            for (int j=0; j < A.rows(); j++) {
                double rowSum = Acache[i].viewRow(j).zSum();
                if (rowSum != 0) {
                    tmpIndices.add(j);
                }
            }

            AcacheIndices[i] = tmpIndices.stream().mapToInt(p -> p).toArray();
            Acache[i] = Acache[i].viewSelection(AcacheIndices[i], null);
            Bcache[i] = B.viewSelection(index.get(i), null);
        }

        for (int i = 0; i < index.size(); i++) {
            DoubleMatrix2D Ctmp = DoubleFactory2D.dense.make(Acache[i].rows(), 1);

            for (int col = 0; col < B.columns(); col++) {

                Acache[i].zMult(Bcache[i].viewSelection(null, new int[]{col}), Ctmp);
                for (int j = 0; j < Ctmp.rows(); j++) {
                    int idx = AcacheIndices[i][j];
                    C.setQuick(idx, col, C.getQuick(idx, col) + Ctmp.get(j, 0));
                }
                pb.step();
            }

        }
        pb.close();

        return C;
    }


    /**
     * Shorthand for A * B
     *
     * @param A A matrix
     * @param B A matrix
     * @return
     */
    private static DoubleMatrix2D mult(DoubleMatrix2D A, DoubleMatrix2D B) {
        //DoubleMatrix2D C = DoubleFactory2D.dense.make(A.rows(), B.columns());
        //A.zMult(B, C);
        // return C;
        return DenseDoubleAlgebra.DEFAULT.mult(A, B);
    }

    private static double mult(DoubleMatrix1D A, DoubleMatrix1D B) {
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
        return DoubleFactory2D.sparse.diagonal(A);
    }

    /**
     * Shorthand for diag(A)
     *
     * @param A
     * @return
     */
    private static DoubleMatrix2D diag(DoubleMatrix1D A) {
        return DoubleFactory2D.sparse.diagonal(A);
    }

    private static DoubleMatrix2D prodDiag(DoubleMatrix2D A, DoubleMatrix1D diag) {

        DoubleMatrix2D results = A.copy();
        for (int j = 0; j < A.columns(); j++) {
            for (int i = 0; i < A.rows(); i++) {
                results.setQuick(i, j, A.getQuick(i, j) * diag.get(i));
            }
        }

        return (results);
    }

    private static DoubleMatrix2D divDiag(DoubleMatrix2D A, DoubleMatrix1D diag) {

        DoubleMatrix2D results = A.copy();
        for (int j = 0; j < A.columns(); j++) {
            for (int i = 0; i < A.rows(); i++) {
                results.setQuick(i, j, A.getQuick(i, j) / diag.get(i));
            }
        }

        return (transpose(results));
    }

    private static DoubleMatrix1D prodDiagT(DoubleMatrix1D A, DoubleMatrix1D diag) {
        DoubleMatrix1D results = DoubleFactory1D.dense.make((int) A.size());
        for (int i = 0; i < A.size(); i++) {
            results.setQuick(i, A.getQuick(i) * diag.get(i));
        }

        return (results);
    }


    private static DoubleMatrix2D prodDiag(DoubleMatrix1D A, DoubleMatrix1D diag) {
        DoubleMatrix2D results = DoubleFactory2D.dense.make((int) A.size(), 1);
        for (int i = 0; i < A.size(); i++) {
            results.setQuick(i, 0, A.getQuick(i) * diag.get(i));
        }

        return (results);
    }



}
