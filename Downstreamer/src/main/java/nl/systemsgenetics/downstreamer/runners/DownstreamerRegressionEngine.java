package nl.systemsgenetics.downstreamer.runners;

import cern.colt.function.tdouble.DoubleDoubleFunction;
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
import org.apache.commons.math3.util.FastMath;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;
import org.jblas.DoubleMatrix;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class DownstreamerRegressionEngine {

    private static final Logger LOGGER = LogManager.getLogger(DownstreamerRegressionEngine.class);

    // Vector functions
    public static final DoubleDoubleFunction minus = (a, b) -> a - b;
    public static final DoubleDoubleFunction plus = (a, b) -> a + b;

    public static void run(OptionsModeRegress options) throws Exception {

        LOGGER.warn("Assumes input data are in the same order unless -ro is specified. Order in -ro should match with -u and -l");
        logInfoMem("Loading datasets");

        DoubleMatrixDataset<String, String> X = DoubleMatrixDataset.loadDoubleData(options.getExplanatoryVariables().getPath());
        logInfoMem("Dim X:" + X.rows() + "x" + X.columns());

        DoubleMatrixDataset<String, String> Y = DoubleMatrixDataset.loadDoubleData(options.getResponseVariable().getPath());
        logInfoMem("Dim Y:" + Y.rows() + "x" + Y.columns());

        // Keep track of overlapping rows
        HashSet<String> overlappingRows = new HashSet<>(X.getRowObjects());
        overlappingRows.retainAll(Y.getRowObjects());

        // Covariates
        DoubleMatrixDataset<String, String> C = null;
        if (options.hasCovariates()) {
            C = DoubleMatrixDataset.loadDoubleData(options.getCovariates().getPath());
            logInfoMem("Dim C:" + C.rows() + "x" + C.columns());

            overlappingRows.retainAll(C.getRowObjects());
        }

        DoubleMatrixDataset<String, String> Sigma = null;
        DoubleMatrixDataset<String, String> U = null;
        DoubleMatrixDataset<String, String> L = null;

        if (options.hasSigma()) {
            throw new Exception("Not yet supported");
        } else {
            U = DoubleMatrixDataset.loadDoubleData(options.getEigenvectors().getPath());
            logInfoMem("Dim U:" + U.rows() + "x" + U.columns());

            L = DoubleMatrixDataset.loadDoubleData(options.getEigenvalues().getPath());
            logInfoMem("Dim L:" + L.rows() + "x" + L.columns());

            overlappingRows.retainAll(U.getRowObjects());
        }

        if (options.hasRowIncludeFilter()) {
            logInfoMem("Filtering rows to rows provided in -ro");
            overlappingRows.retainAll(IoUtils.readMatrixAnnotations(options.getRowIncludeFilter()));
        }

        if (options.hasColumnIncludeFilter()) {
            logInfoMem("Filtering cols in X to those provided in -co");
            HashSet<String> colsToSelect = new HashSet<>(IoUtils.readMatrixAnnotations(options.getColumnIncludeFilter()));
            X = X.viewColSelection(colsToSelect);
        }

        // Scanning for NA's
        logInfoMem("Screening for NA values. These are removed row wise. If you want to remove them column wise, please do so first.");

        Set<String> rowsWithNA = checkNaRowWise(Y);
        rowsWithNA.addAll(checkNaRowWise(Y));

        if (C != null) {
            rowsWithNA.addAll(checkNaRowWise(C));
        }

        if (U != null) {
            rowsWithNA.addAll(checkNaRowWise(U));
        }

        if (Sigma != null) {
            rowsWithNA.addAll(checkNaRowWise(Sigma));
        }

        overlappingRows.removeAll(rowsWithNA);
        logInfoMem("Removed " + rowsWithNA.size() + " rows with NA values.");

        // Convert to list to ensure consistent order
        //final List<String> finalRowSelection = new ArrayList<>(overlappingRows);

        List<String> finalRowSelection = Y.getRowObjects();
        finalRowSelection.remove(0);

        // Subset on overlapping rows and make sure everything is in the same order
        X = X.viewRowSelection(finalRowSelection);
        Y = Y.viewRowSelection(finalRowSelection);

        if (C != null) {
            C = C.viewRowSelection(finalRowSelection);
        }

        if (U != null) {
            U = U.viewRowSelection(finalRowSelection);
        }

        if (Sigma != null) {
            Sigma = Sigma.viewSelection(finalRowSelection, finalRowSelection);
        }

        logInfoMem("Done loading datasets");
        logInfoMem("Maintaining " + overlappingRows.size() + " overlapping rows");


        // Inverse normal transform
        if (options.isInverseNormalY()) {
            logInfoMem("Applying INT to Y");
            Y.createColumnForceNormalInplace();
        }

        if (options.isInverseNormalX()) {
            logInfoMem("Applying INT to X");
            X.createColumnForceNormalInplace();
        }

        if (options.isInverseNormalC() && C != null) {
            logInfoMem("Applying INT to C");
            C.createColumnForceNormalInplace();
        }

        // Center and scale.
        if (options.centerAndScale()) {
            centerAndScale(Y);
            centerAndScale(X);
            if (C != null) {
                centerAndScale(C);
            }
        }

        if (options.regressCovariates()) {
            logInfoMem("Regressing out effects of covariates before main regression using OLS model.");
            inplaceDetermineOlsRegressionResiduals(Y, C);
            // Set covariates to null so they are not used for the rest of the regression
            C = null;
        }

        // At this point it is assumed all matching of genes has been done. There should be no missing genes in
        // the genes file (this is checked and error is thrown). If there are, users can use -ro
        List<int[]> blockDiagonalIndices = null;

        if (options.hasGenes()) {
            blockDiagonalIndices = createBlockDiagonalIndexFromGenes(options.getGenes(), X.getRowObjects());
        }

        // Depending on input, first decompse Sigma, or run with precomputed eigen decomp.
        LinearRegressionResult output;
        if (options.hasSigma()) {
            throw new Exception("Not supported atm. Please pre-calculate the eigen decomposition");
        } else {
            output = performDownstreamerRegression(X, Y, C, U, L,
                    blockDiagonalIndices,
                    options);
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
                                                                       OptionsModeRegress options) {

        DoubleMatrix1D varPerEigen = cumulativeSum(L.viewCol(0), true);
        List<String> colsToMaintain = new ArrayList<>();

        for (int i = 0; i < varPerEigen.size(); i++) {
            if (options.getPercentageOfVariance() > varPerEigen.get(i)) {
                colsToMaintain.add(U.getColObjects().get(i));
            }
        }

        int numVectorsToMaintain = colsToMaintain.size();
        logInfoMem("Maintaining " + numVectorsToMaintain + " eigenvectors explaining " + varPerEigen.get(numVectorsToMaintain - 1) + " % of the variance.");

        // Subset U and L on the number of eigenvectors to maintain
        U = U.viewColSelection(colsToMaintain);
        L = L.viewRowSelection(colsToMaintain);

        // Pre-transpose U as this can be re-used as well
        DoubleMatrix2D UHatT = transpose(U.getMatrix());
        //UHatT = DoubleFactory2D.sparse.make(UHatT.toArray());

        // Calculate inverse of eigenvalues
        DoubleMatrix1D LHatInv = L.viewCol(0).copy();
        LHatInv.assign(DoubleFunctions.inv);

        // Pre-compute Y^ as this can be re-used
        DoubleMatrix2D YHat = mult(UHatT, Y.getMatrix());

        logInfoMem("Starting regression for " + X.columns() + " pathways.");

        // Determine the degrees of freedom. -1 for main_effect
        int degreesOfFreedom = numVectorsToMaintain - 1;

        // Set the names of the predictors in the correct order
        // and determine the degrees of freedom
        List<String> predictorNames = new ArrayList<>();
        if (options.fitIntercept()) {
            predictorNames.add("intercept");
            degreesOfFreedom = degreesOfFreedom - 1;
        }
        predictorNames.add("main_effect");
        if (C != null) {
            predictorNames.addAll(C.getColObjects());
            degreesOfFreedom = degreesOfFreedom - C.columns();
        }

        // Pre calculate predictors and covariates in eigenvectors space.
        DoubleMatrix2D XhatCache = null;
        DoubleMatrix2D ChatCache = null;

        if (blockDiagonalIndices != null) {
            logInfoMem("Precomputing XHat and CHat");
            XhatCache = multBlockDiagonalSubsetPerCol(UHatT, X.getMatrix(), blockDiagonalIndices, options.useJblas());

            if (C != null) {
                ChatCache = multBlockDiagonalSubsetPerCol(UHatT, C.getMatrix(), blockDiagonalIndices, options.useJblas());
            }
            logInfoMem("Done");
        } else {
            if (C != null) {
                ChatCache = mult(UHatT, C.getMatrix());
            }
        }

        // Object to save output
        LinearRegressionResult result = new LinearRegressionResult(X.getColObjects(), predictorNames, degreesOfFreedom);
        ProgressBar pb = new ProgressBar("Linear regressions", X.columns());

        // Calculate beta's and SE for each pathway
        for (int curPathway = 0; curPathway < X.columns(); curPathway++) {

            // Get the current pathway, if not pre-multiplied by eigenvectors, do so now.
            DoubleMatrix2D XCur;
            if (XhatCache == null) {
                XCur = mult(UHatT, X.viewColAsMmatrix(curPathway));
            } else {
                XCur = XhatCache.viewSelection(null, new int[]{curPathway});
            }

            // If covariates are provided add them
            if (ChatCache != null) {
                XCur = DoubleFactory2D.dense.appendColumns(XCur, ChatCache);
            }

            // If not pre-computing Xhat, center and scale here
            if (options.centerAndScale() && XhatCache == null) {
                centerAndScale(XCur);
            }

            double[] curRes = downstreamerRegressionPrecomp(XCur,
                    YHat,
                    LHatInv,
                    options.fitIntercept());

            result.appendBetas(curPathway, ArrayUtils.subarray(curRes, 0, curRes.length / 2));
            result.appendSe(curPathway, ArrayUtils.subarray(curRes, curRes.length / 2, curRes.length));

            pb.step();
        }
        pb.close();

        logInfoMem("Done with regression for " + X.columns() + " pathways.");

        return result;
    }


    /**
     * Generic implementation of the downstreamer regression model where:
     * XHat = t(UHat) * X
     * YHat = t(UHat) * Y
     * LHatInv = 1/diag(LHat)
     * have been pre-computed to save time as these remain constant for all pathways in X.
     *
     * @param XHat         Predictors, e x p matrix XHat = t(UHat) * X
     * @param YHat         Response, e x 1 matrix t(UHat) * Y
     * @param LHatInv      Inverse eigenvalues, e x 1 matrix, 1/diag(LHat)
     * @param fitIntercept Should an intercept term be fitted
     * @return double[] with beta and standard erros in the form  [beta_1, beta_2 ... beta_x, se_1, se_2 ... se_x]
     */
    public static double[] downstreamerRegressionPrecomp(DoubleMatrix2D XHat, DoubleMatrix2D YHat, DoubleMatrix1D LHatInv, boolean fitIntercept) {

        // TODO: DoubleMatrix2D
        double[] results = new double[XHat.columns() * 2];

        if (fitIntercept) {
            results = new double[(XHat.columns() * 2) + 2];
            XHat = DoubleFactory2D.dense.appendColumns(DoubleFactory2D.dense.make(XHat.rows(), 1, 1), XHat);
        }

        //-----------------------------------------------------------------------------------------
        // beta
        //-----------------------------------------------------------------------------------------
        //   crossprod = t(A) * B > this is the R code used for testing
        //  a    <- crossprod(Y.hat, L.hat.inv) %*% X.hat
        //  b    <- solve(crossprod(X.hat, L.hat.inv) %*% X.hat)
        //  beta <- a %*% b

        //DoubleMatrix2D a = mult(tmult(YHat, LHatInv), XHat);
        //DoubleMatrix2D b = inverse(mult(tmult(XHat, LHatInv), XHat));

        // Alternative for YHat %*% LHatInv, bit quicker
        DoubleMatrix2D a = mult(transpose(multDiag(YHat, LHatInv)), XHat);
        DoubleMatrix2D b = inverse(mult(transpose(multDiag(XHat, LHatInv)), XHat));
        DoubleMatrix2D beta = mult(a, b);

        //-----------------------------------------------------------------------------------------
        // Residuals
        //-----------------------------------------------------------------------------------------
        // Predicted Y
        DoubleMatrix1D predicted = multT(XHat, beta).viewColumn(0);
        // Original Y
        DoubleMatrix1D residuals = YHat.viewColumn(0).copy();
        // Subtract predicted from original
        residuals.assign(predicted, minus);

        // Residual sum of squares
        //double rss = mult(tmult(residuals, LHatInv), residuals.reshape((int) residuals.size(), 1)).get(0, 0);

        // Alternative for residuals %*% LHatInv, bit quicker,
        double rss = mult(tMultDiag(residuals, LHatInv), residuals);

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

    /**
     * Checks for NA values, returns an index of which rows have NA.
     */
    public static Set<String> checkNaRowWise(DoubleMatrixDataset<String, String> input) {

        Set<String> rowsWithNa = new HashSet<>();

        for (int i = 0; i < input.rows(); i++) {
            for (int j = 0; j < input.columns(); j++) {
                if (Double.isNaN(input.getElementQuick(i, j))) {
                    rowsWithNa.add(input.getRowObjects().get(i));
                    break;
                }
            }
        }
        return rowsWithNa;
    }


    /**
     * Creates a list of row indices that are on a chromsome arm, and are thus block diagonal.
     * This is assumed to be true and is not checked.
     *
     * @param geneFile Ensembl file specified by --genes
     * @param rowOrder A list with gene ids in the order of the data.
     * @return
     * @throws IOException
     */
    private static List<int[]> createBlockDiagonalIndexFromGenes(File geneFile, List<String> rowOrder) throws IOException {

        // Determine the indices in the matrix that lie together on a chromosome arm
        List<int[]> blockDiagonalIndices = new ArrayList<>();
        Map<String, List<Integer>> blockDiagonalIndicesTmp = new HashMap<>();
        Map<String, List<Gene>> genes = IoUtils.readGenesAsChrArmMap(geneFile);

        int totalCount = 0;

        for (String key : genes.keySet()) {
            for (Gene curGene : genes.get(key)) {
                int idx = rowOrder.indexOf(curGene.getGene());
                if (idx > -1) {
                    if (!blockDiagonalIndicesTmp.containsKey(key)) {
                        blockDiagonalIndicesTmp.put(key, new ArrayList<>());
                    }
                    blockDiagonalIndicesTmp.get(key).add(idx);
                    totalCount++;
                }
            }
        }

        if (rowOrder.size() > totalCount) {
            throw new IllegalArgumentException("There were genes on rows not found in --genes file. " +
                    "Please make sure all genes are available in --genes or use -ro to make select a subset.");
        }

        // Convert to int[] so it's easier to use with DoubleMatrix2D;
        for (List<Integer> curList : blockDiagonalIndicesTmp.values()) {
            blockDiagonalIndices.add(curList.stream().mapToInt(i -> i).toArray());
        }

        return blockDiagonalIndices;
    }


    // Arithmatic functions to make code more readable //

    private static void centerAndScale(DoubleMatrixDataset matrix) {
        centerAndScale(matrix.getMatrix());
    }


    /**
     * Set mean to 0 and sd to 1 for each column in the data.
     * Lifted from DoubleMatrixDataset.normalizeColumns()
     * Refactored to center and scale to avoid confusion with DoubleMatrix2D.normalize()
     * which normalizes the data so their sum = 1
     *
     * @param matrix
     */
    private static void centerAndScale(DoubleMatrix2D matrix) {

        final int rows = matrix.rows();

        for (int c = 0; c < matrix.columns(); ++c) {

            DoubleMatrix1D col = matrix.viewColumn(c);

            double colSum = 0;

            for (int e = 0; e < rows; ++e) {
                colSum += col.getQuick(e);
            }

            final double mean = colSum / rows;

            double varSum = 0;
            for (int e = 0; e < rows; ++e) {
                varSum += (col.getQuick(e) - mean) * (col.getQuick(e) - mean);
            }

            double sd = FastMath.sqrt(varSum / (rows - 1));

            for (int e = 0; e < rows; ++e) {
                col.setQuick(e, (col.getQuick(e) - mean) / sd);
            }

        }

    }


    /**
     * Calculate the cumlative sum of the elements in A.
     *
     * @param a            Vector to calculate cumulative sum for
     * @param asPercentage Should results be returned as a proportion
     * @return
     */
    private static DoubleMatrix1D cumulativeSum(DoubleMatrix1D a, boolean asPercentage) {

        double totalSum = 1;
        if (asPercentage) {
            totalSum = a.zSum();
        }
        DoubleMatrix1D result = a.like();
        double runningSum = 0;
        for (int i = 0; i < a.size(); i++) {
            runningSum = runningSum + a.get(i);
            result.set(i, runningSum / totalSum);
        }

        return result;
    }


    /**
     * Calculate A %*% B for a matrix A whose structure has some form of block diagonality. I.e. there are rows
     * and columns that are exactly zero. @param index gives the indices to multiply together.
     * JBlas is much quicker on larger matrices but has memory overhead compared to colt matrix mult.
     *
     * @param A
     * @param B
     * @param index
     * @param useJblas
     * @return
     */
    private static DoubleMatrix2D multBlockDiagonalSubsetPerCol(DoubleMatrix2D A, DoubleMatrix2D B, List<int[]> index, boolean useJblas) {
        if (useJblas) {
            return multBlockDiagonalSubsetPerColJblas(A, B, index);
        } else {
            return multBlockDiagonalSubsetPerColPcolt(A, B, index);
        }
    }

    /**
     * Calculate A %*% B for a matrix A whose structure has some form of block diagonality. I.e. there are rows
     * and columns that are exactly zero. @param index gives the indices to multiply together.
     * Uses Jblas for matrix multiplications. This does give memory overhead but is much quciker.
     *
     * @param Ac
     * @param Bc
     * @param index
     * @return
     */
    private static DoubleMatrix2D multBlockDiagonalSubsetPerColJblas(DoubleMatrix2D Ac, DoubleMatrix2D Bc, List<int[]> index) {

        DoubleMatrix2D C = DoubleFactory2D.dense.make(Ac.rows(), Bc.columns());
        ProgressBar pb = new ProgressBar("Multiplying per block", Bc.columns() * index.size());

        DoubleMatrix[] Acache = new DoubleMatrix[index.size()];
        DoubleMatrix[] Bcache = new DoubleMatrix[index.size()];

        int[][] AcacheIndices = new int[index.size()][];

        // Compile the subset caches. Precomputing saves overhead, but is more memory intense.
        for (int i = 0; i < index.size(); i++) {
            Acache[i] = new DoubleMatrix(Ac.viewSelection(null, index.get(i)).toArray());
            DoubleMatrix tmp = Acache[i].rowSums();
            ArrayList<Integer> tmpIndices = new ArrayList<>();
            for (int j = 0; j < tmp.rows; j++) {
                if (tmp.get(j, 0) != 0) {
                    tmpIndices.add(j);
                }
            }

            AcacheIndices[i] = tmpIndices.stream().mapToInt(p -> p).toArray();
            Acache[i] = Acache[i].getRows(AcacheIndices[i]);
            Bcache[i] = new DoubleMatrix(Bc.viewSelection(index.get(i), null).toArray());
        }

        // Perform per block multiplication
        for (int i = 0; i < index.size(); i++) {
            // Re-cyle temorary results matrix
            DoubleMatrix Ctmp = new DoubleMatrix(Acache[i].rows, 1);

            for (int col = 0; col < Bc.columns(); col++) {
                // Jblas multiplication
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

    /**
     * Calculate A %*% B for a matrix A whose structure has some form of block diagonality. I.e. there are rows
     * and columns that are exactly zero. @param index gives the indices to multiply together.
     *
     * @param Ac
     * @param Bc
     * @param index
     * @return
     */
    private static DoubleMatrix2D multBlockDiagonalSubsetPerColPcolt(DoubleMatrix2D A, DoubleMatrix2D B, List<int[]> index) {

        DoubleMatrix2D C = DoubleFactory2D.dense.make(A.rows(), B.columns());

        ProgressBar pb = new ProgressBar("Multiplying per block", B.columns() * index.size());

        DoubleMatrix2D[] Acache = new DoubleMatrix2D[index.size()];
        DoubleMatrix2D[] Bcache = new DoubleMatrix2D[index.size()];

        int[][] AcacheIndices = new int[index.size()][];

        // Compile the subset caches. Precomputing saves overhead. Views are returned so is memory efficient.
        for (int i = 0; i < index.size(); i++) {
            Acache[i] = A.viewSelection(null, index.get(i));
            ArrayList<Integer> tmpIndices = new ArrayList<>();

            for (int j = 0; j < A.rows(); j++) {
                double rowSum = Acache[i].viewRow(j).zSum();
                if (rowSum != 0) {
                    tmpIndices.add(j);
                }
            }

            AcacheIndices[i] = tmpIndices.stream().mapToInt(p -> p).toArray();
            Acache[i] = Acache[i].viewSelection(AcacheIndices[i], null);
            Bcache[i] = B.viewSelection(index.get(i), null);
        }

        // Perform per block multiplication
        for (int i = 0; i < index.size(); i++) {
            // Re-cyle temorary results matrix
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
     * Fit an OLS model to each column in Y seperately and overwrite the values of Y by the residuals of the regression.
     * @param Y
     * @param X
     */
    public static void inplaceDetermineOlsRegressionResiduals(DoubleMatrixDataset Y, DoubleMatrixDataset X) {

        // Build design matrix with intercept
        DoubleMatrix2D design = DoubleFactory2D.dense.appendColumns(DoubleFactory2D.dense.make(X.rows(), 1, 1), X.getMatrix());

        // b = (t(X) X)^-1 t(X) y
        // Pre compute this component as thsi is constant for the design matrix
        DoubleMatrix2D XtXi = inverse(tMult(design, design));

        // Loop over columns in y
        for (int col=0; col < Y.columns(); col++) {

            // Beta
            DoubleMatrix2D y = Y.viewCol(col).reshape(Y.rows(), 1);
            DoubleMatrix2D betas = mult(mult(XtXi, design.viewDice()), y);

            // Predicted Y
            DoubleMatrix2D predicted = mult(betas.viewDice(), design.viewDice());

            // Subtract predicted from original to get residuals
            Y.viewCol(col).assign(predicted.viewRow(0), minus);
        }
    }


    /**
     * Shorthand for A * t(B)
     *
     * @param A A matrix
     * @param B A matrix
     * @return
     */
    private static DoubleMatrix2D multT(DoubleMatrix2D A, DoubleMatrix2D B) {
        return mult(A, B.viewDice());
    }

    /**
     * Shorthand for t(A) * B
     *
     * @param A A matrix
     * @param B A matrix
     * @return
     */
    private static DoubleMatrix2D tMult(DoubleMatrix2D A, DoubleMatrix2D B) {
        return mult(A.viewDice(), B);
    }

    /**
     * Shorthand for t(A) * B
     *
     * @param A A column vector
     * @param B A matrix
     * @return
     */
    private static DoubleMatrix2D tMult(DoubleMatrix1D a, DoubleMatrix2D B) {
        return mult(a.reshape((int) a.size(), 1).viewDice(), B);
    }

    /**
     * Shorthand for A * B but save the result in C
     *
     * @param A A matrix
     * @param B A matrix
     * @return
     */
    private static DoubleMatrix2D multI(DoubleMatrix2D A, DoubleMatrix2D B, DoubleMatrix2D C) {
        return A.zMult(B, C);
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

    private static double mult(DoubleMatrix1D a, DoubleMatrix1D b) {
        return DenseDoubleAlgebra.DEFAULT.mult(a, b);
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
        return A.viewDice();
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
    private static DoubleMatrix2D diag(DoubleMatrix1D a) {
        return DoubleFactory2D.sparse.diagonal(a);
    }

    /**
     * Calculate the pairwise product between a matrix and a "diagonal" matrix column wise. The diagonal is represented
     * as a vector.
     *
     * @param A
     * @param diag
     * @return
     */
    private static DoubleMatrix2D multDiag(DoubleMatrix2D A, DoubleMatrix1D diag) {
        DoubleMatrix2D results = A.copy();
        for (int j = 0; j < A.columns(); j++) {
            for (int i = 0; i < A.rows(); i++) {
                results.setQuick(i, j, A.getQuick(i, j) * diag.get(i));
            }
        }
        return (results);
    }

    /**
     * Calculate the pairwise product between a vector and a diagonal matrix column wise. The diagonal is represented
     * as a vector.
     *
     * @param A
     * @param diag
     * @return
     */
    private static DoubleMatrix1D tMultDiag(DoubleMatrix1D a, DoubleMatrix1D diag) {
        DoubleMatrix1D results = DoubleFactory1D.dense.make((int) a.size());
        for (int i = 0; i < a.size(); i++) {
            results.setQuick(i, a.getQuick(i) * diag.get(i));
        }

        return (results);
    }

    private static void logInfoMem(String msg) {
        Downstreamer.logInfoMem(msg, LOGGER);
    }

}
