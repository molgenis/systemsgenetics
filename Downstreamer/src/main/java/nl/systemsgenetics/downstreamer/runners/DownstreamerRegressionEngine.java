package nl.systemsgenetics.downstreamer.runners;

import cern.colt.function.tdouble.DoubleDoubleFunction;
import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import cern.colt.matrix.tdouble.algo.decomposition.DenseDoubleEigenvalueDecomposition;
import cern.jet.math.tdouble.DoubleFunctions;
import gnu.trove.list.array.TIntArrayList;
import me.tongfei.progressbar.ProgressBar;
import nl.systemsgenetics.downstreamer.Downstreamer;
import nl.systemsgenetics.downstreamer.gene.Gene;
import nl.systemsgenetics.downstreamer.gene.IndexedDouble;
import nl.systemsgenetics.downstreamer.io.BlockDiagonalDoubleMatrixProvider;
import nl.systemsgenetics.downstreamer.io.DoubleMatrixDatasetBlockDiagonalProvider;
import nl.systemsgenetics.downstreamer.io.IoUtils;
import nl.systemsgenetics.downstreamer.runners.options.OptionsModeRegress;
import nl.systemsgenetics.downstreamer.summarystatistic.LinearRegressionResult;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.NotImplementedException;
import org.apache.commons.math3.util.FastMath;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;
import org.jblas.DoubleMatrix;
import org.jblas.Eigen;
import smile.stat.Hypothesis;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class DownstreamerRegressionEngine {

	private static final Logger LOGGER = LogManager.getLogger(DownstreamerRegressionEngine.class);

	// Vector functions
	public static final DoubleDoubleFunction minus = (a, b) -> a - b;
	public static final DoubleDoubleFunction plus = (a, b) -> a + b;

	public static void run(OptionsModeRegress options) throws Exception {

		logInfoMem("Loading datasets");

		DoubleMatrixDataset<String, String> X = DoubleMatrixDataset.loadDoubleData(options.getExplanatoryVariables().getPath());
		logInfoMem("Dim X:" + X.rows() + "x" + X.columns());

		DoubleMatrixDataset<String, String> Y = DoubleMatrixDataset.loadDoubleData(options.getResponseVariable().getPath());
		logInfoMem("Dim Y:" + Y.rows() + "x" + Y.columns());

		// Keep track of overlapping rows
		LinkedHashSet<String> overlappingRows = new LinkedHashSet<>(X.getRowObjects());
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
            //throw new Exception("Eigen decomp of sigma is not yet supported here, please pre-compute and use -u -l");
            Sigma = DoubleMatrixDataset.loadDoubleData(options.getSigma().getPath());
            logInfoMem("Dim S:" + Sigma.rows() + "x" + Sigma.columns());

            overlappingRows.retainAll(Sigma.getRowObjects());
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

		// Convert to list to ensure consistent order, altough might be redundant as linkedHasSet should respect order
		final List<String> finalRowSelection = new ArrayList<>(overlappingRows);

		// Subset on overlapping rows and make sure everything is in the same order
		X = X.viewRowSelection(finalRowSelection);
		Y = Y.viewRowSelection(finalRowSelection);

		if (C != null) {
			C = C.viewRowSelection(finalRowSelection);
		}

        if (U != null) {
            U = U.viewRowSelection(finalRowSelection);

            if (U.columns() > 1.1 * U.rows()) {
                LOGGER.warn("!!!MORE THAN 10% DIFFERENCE BETWEEN THE NUMBER OF GENES AND NUMBER OF EIGENVECTORS!!!");
                LOGGER.warn("This will lead to an overestimation of the degrees of freedom and produce invalid standard errors");
                LOGGER.warn("We strongly advise to calculate the eigenvectors only on genes that overlap!");
            }
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

		// At this point it is assumed all matching of genes has been done. There should be no missing genes in
		// the genes file (this is checked and error is thrown). If there are, users can use -ro
		List<int[]> blockDiagonalIndices = null;

        if (options.hasGenes()) {
            blockDiagonalIndices = createBlockDiagonalIndexFromGenes(options.getGenes(), finalRowSelection);

            if (options.isDebugMode()) {
                IoUtils.writeBlockDiagonalIndices(blockDiagonalIndices, finalRowSelection, options.getDebugFolder() + "/block_diagonal_indices.txt");
            }
        }

        // Depending on input, first decompse Sigma, or run with precomputed eigen decomp.
        List<LinearRegressionResult> output;

        // TODO: in future make sure this doesn't read the full matrix
        BlockDiagonalDoubleMatrixProvider sigmaProvider = new DoubleMatrixDatasetBlockDiagonalProvider(Sigma);
        DoubleMatrixDataset<String, String>[] eigen = null;

        if (Sigma != null && blockDiagonalIndices != null) {
            eigen = blockDiagonalEigenDecomposition(sigmaProvider, blockDiagonalIndices, options.useJblas());
        } else if (Sigma != null) {
            eigen = eigenDecomposition(Sigma, options.useJblas());
        }

        if (eigen != null) {
            L = eigen[0];
            U = eigen[1];

            U = U.viewRowSelection(finalRowSelection);

            if (options.isDebugMode()) {
                U.save(options.getDebugFolder() + "/genecor_eigenvectors.txt");
                L.save(options.getDebugFolder() + "/genecor_eigenvalues.txt");

				BufferedWriter writer = new BufferedWriter(new FileWriter(options.getDebugFolder() + "/final_row_order.txt"));
				for (String cur: finalRowSelection) {
					writer.write(cur);
					writer.newLine();
				}
				writer.flush();
				writer.close();
            }
        }

        output = performDownstreamerRegression(X, Y, C, U, L,
                blockDiagonalIndices,
                options.getPercentageOfVariance(),
                options.fitIntercept(),
                options.regressCovariates(),
                options.useJblas());

        // Save linear regression results
        for (LinearRegressionResult curRes : output) {
            curRes.save(options.getOutputBasePath() + "_" + curRes.getName(), false);
        }
    }

	/**
	 * Runs DS regression if eigendecompostion has been pre-computed.These are
	 * given by U and L. Expects that these still need to be filtered. If they
	 * are pre-selected make sure to set percentageOfVariance to 1 or call other
	 * implementation below.
	 *
	 * @param X pathway scores
	 * @param Y trait gene scores
	 * @param C covariates on columns genes on rows
	 * @param U eigenvectors that describe gene-gene correlations
	 * @param L single column with on each row the eigenvalue of a component
	 * @param blockDiagonalIndices
	 * @param percentageOfVariance determines the number of eigenvectors to use.
	 * Set to 1 if you want to use all
	 * @param fitIntercept should an intercept term be fitted
	 * @param regressCovariates should covariates be regressed first instead of
	 * fit in the model with X
	 * @param useJblas use Jblas non-java matrix multiplications instead of
	 * Pcolt. Faster but may give system specific issues
	 * @return Per trait (column in Y) a LinearRegressionResult.
	 */
	public static List<LinearRegressionResult> performDownstreamerRegression(DoubleMatrixDataset<String, String> X,
			DoubleMatrixDataset<String, String> Y,
			DoubleMatrixDataset<String, String> C,
			DoubleMatrixDataset<String, String> U,
			DoubleMatrixDataset<String, String> L,
			List<int[]> blockDiagonalIndices,
			double percentageOfVariance,
			boolean fitIntercept,
			boolean regressCovariates,
			boolean useJblas) {

		DoubleMatrix1D varPerEigen = cumulativeSum(L.viewCol(0), true);
		List<String> colsToMaintain = new ArrayList<>();

		for (int i = 0; i < varPerEigen.size(); i++) {
			if (percentageOfVariance > varPerEigen.get(i)) {
				colsToMaintain.add(U.getColObjects().get(i));
			}
		}

		int numVectorsToMaintain = colsToMaintain.size();
		logInfoMem("Maintaining " + numVectorsToMaintain + " eigenvectors explaining " + varPerEigen.get(numVectorsToMaintain - 1) + " % of the variance.");

		// Subset U and L on the number of eigenvectors to maintain
		U = U.viewColSelection(colsToMaintain);
		L = L.viewRowSelection(colsToMaintain);

		return performDownstreamerRegression(X, Y, C, U, L, blockDiagonalIndices, fitIntercept, regressCovariates, useJblas);
	}

	/**
	 * Runs DS regression if U and L have been pre-filtered on only the
	 * eigenvectors to use. See above for parameter descriptions as they are
	 * identical.
	 */
	public static List<LinearRegressionResult> performDownstreamerRegression(DoubleMatrixDataset<String, String> X,
			DoubleMatrixDataset<String, String> Y,
			DoubleMatrixDataset<String, String> C,
			DoubleMatrixDataset<String, String> U,
			DoubleMatrixDataset<String, String> L,
			List<int[]> blockDiagonalIndices,
			boolean fitIntercept,
			boolean regressCovariates,
			boolean useJblas) {

		// Pre-transpose U as this can be re-used as well
		DoubleMatrix2D UHatT = transpose(U.getMatrix());
		//UHatT = DoubleFactory2D.sparse.make(UHatT.toArray());

		// Calculate inverse of eigenvalues
		DoubleMatrix1D LHatInv = L.viewCol(0).copy();
		LHatInv.assign(DoubleFunctions.inv);

		logInfoMem("Starting regression for " + X.columns() + " pathways.");

        // Determine the degrees of freedom. -1 for main_effect
        int degreesOfFreedom = UHatT.rows() - 1;

		// Set the names of the predictors in the correct order
		// and determine the degrees of freedom
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

		// Pre calculate predictors and covariates in eigenvectors space.
		DoubleMatrix2D XhatCache = null;
		DoubleMatrix2D ChatCache = null;

		// Pre-calculate XHat and Xhat in a block diagonal way, or using full matrix mult.
		if (blockDiagonalIndices != null) {
			logInfoMem("Precomputing XHat and CHat");
			XhatCache = blockDiagonalMult(UHatT, X.getMatrix(), blockDiagonalIndices, useJblas);

            if (C != null) {
                ChatCache = blockDiagonalMult(UHatT, C.getMatrix(), blockDiagonalIndices, useJblas);
            }
            logInfoMem("Done");
        } else {
            if (C != null) {
                ChatCache = mult(UHatT, C.getMatrix());
            }
        }

		ProgressBar pb = new ProgressBar("Linear regressions", X.columns() * Y.columns());

		// Store output
		List<LinearRegressionResult> results = new ArrayList<>(Y.columns());

		for (int curY = 0; curY < Y.columns(); curY++) {

			// Keep track of current covariates, if these are used to regress need to keep ChatCache.
			// but curCovariates can be set to null
			DoubleMatrix2D curCovariates = ChatCache;
			List<String> curPredictorNames = predictorNames;

			// Pre-compute Y^ as this can be re-used
			DoubleMatrix2D YHat = mult(UHatT, Y.viewColAsMmatrix(curY));

			// Instead of including covariates in the model, first regress their effects
			if (regressCovariates) {
				inplaceDownstreamerRegressionResidualsPrecomp(YHat, ChatCache, LHatInv);

				// Set covariates to null so they are not used for the rest of the regression
				curCovariates = null;
				curPredictorNames = predictorNames.subList(0, predictorNames.size() - C.columns());
			}

			// Object to save output
			LinearRegressionResult result = new LinearRegressionResult(X.getColObjects(),
					curPredictorNames,
					degreesOfFreedom,
					Y.getColObjects().get(curY));

			// Calculate beta's and SE for each pathway
			for (int curPathway = 0; curPathway < X.columns(); curPathway++) {

				// Get the current pathway, if not pre-multiplied by eigenvectors, do so now.
				DoubleMatrix2D XCur;
				if (XhatCache == null) {
					XCur = mult(UHatT, X.viewColAsMmatrix(curPathway));
				} else {
					XCur = XhatCache.viewSelection(null, new int[]{curPathway});
				}

				// If covariates are provided add them to the design matrix
				if (curCovariates != null) {
					XCur = DoubleFactory2D.dense.appendColumns(XCur, ChatCache);
				}

				double[] curRes = downstreamerRegressionPrecomp(XCur,
						YHat,
						LHatInv,
						fitIntercept);

				result.appendBetas(curPathway, ArrayUtils.subarray(curRes, 0, curRes.length / 2));
				result.appendSe(curPathway, ArrayUtils.subarray(curRes, curRes.length / 2, curRes.length));

				pb.step();
			}

			results.add(result);
		}

		pb.close();
		logInfoMem("Done with regression for " + Y.columns() + " traits and " + X.columns() + " pathways.");

		return results;
	}

	/**
	 * Generic implementation of the downstreamer regression model where: XHat =
	 * t(UHat) * X YHat = t(UHat) * Y LHatInv = 1/diag(LHat) have been
	 * pre-computed to save time as these remain constant for all pathways in X.
	 *
	 * @param XHat Predictors, e x p matrix XHat = t(UHat) * X
	 * @param YHat Response, e x 1 matrix t(UHat) * Y
	 * @param LHatInv Inverse eigenvalues, e x 1 matrix, 1/diag(LHat)
	 * @param fitIntercept Should an intercept term be fitted
	 * @return double[] with beta and standard erros in the form [beta_1, beta_2
	 * ... beta_x, se_1, se_2 ... se_x]
	 */
	private static double[] downstreamerRegressionPrecomp(DoubleMatrix2D XHat, DoubleMatrix2D YHat, DoubleMatrix1D LHatInv, boolean fitIntercept) {

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
	 * Fit a regression model in eigenvector space using precomputed components.
	 * Residuals replace original values in Y Models are fit per col of Y. This
	 * function is usefull for regressing the effects of covariates.
	 *
	 * @param Y Response variables, rows are samples. Should be in eigenvector
	 * space
	 * @param X Predictors to regress, rows are samples. Should be in
	 * eigenvectors space
	 * @param LHatInv Matrix containing the eigenvalues.
	 */
	private static void inplaceDownstreamerRegressionResidualsPrecomp(DoubleMatrix2D Y, DoubleMatrix2D X, DoubleMatrix1D LHatInv) {

		DoubleMatrix2D design = DoubleFactory2D.dense.appendColumns(DoubleFactory2D.dense.make(X.rows(), 1, 1), X);

		// For docs on the stats see comments in downstreamerRegressionPrecomp() above
		DoubleMatrix2D b = inverse(mult(transpose(multDiag(design, LHatInv)), design));

        // Loop over columns in y
        for (int col = 0; col < Y.columns(); col++) {

			DoubleMatrix2D a = mult(transpose(multDiag(Y.viewColumn(col).reshape(Y.rows(), 1), LHatInv)), design);
			DoubleMatrix2D betas = mult(a, b);

			// Predicted Y
			DoubleMatrix1D predicted = multT(design, betas).viewColumn(0);

			// Subtract predicted from original to get residuals
			Y.viewColumn(col).assign(predicted, minus);
		}
	}

	/**
	 * Fit an OLS model to each column in Y seperately and overwrite the values
	 * of Y by the residuals of the regression.
	 *
	 * @param Y Response variables, rows are samples.
	 * @param X Predictors to regress, rows are samples.
	 */
	private static void inplaceDetermineOlsRegressionResiduals(DoubleMatrixDataset Y, DoubleMatrixDataset X) {

		// Build design matrix with intercept
		DoubleMatrix2D design = DoubleFactory2D.dense.appendColumns(DoubleFactory2D.dense.make(X.rows(), 1, 1), X.getMatrix());

		// b = (t(X) X)^-1 t(X) y
		// Pre compute this component as thsi is constant for the design matrix
		DoubleMatrix2D XtXi = inverse(tMult(design, design));

        // Loop over columns in y
        for (int col = 0; col < Y.columns(); col++) {

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
	 * Checks for NA values, returns an index of which rows have NA.
	 */
	private static Set<String> checkNaRowWise(DoubleMatrixDataset<String, String> input) {

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
	 * Creates a list of row indices that are on a chromsome arm, and are thus
	 * block diagonal. This is assumed to be true and is not checked.
	 *
	 * @param geneFile Ensembl file specified by --genes
	 * @param rowOrder A list with gene ids in the order of the data.
	 * @return
	 * @throws IOException
	 */
	public static List<int[]> createBlockDiagonalIndexFromGenes(File geneFile, List<String> rowOrder) throws IOException {
		Map<String, List<Gene>> genes = IoUtils.readGenesAsChrArmMap(geneFile);
		return createBlockDiagonalIndexFromGenes(genes, rowOrder);
	}

	/**
	 * Creates a list of row indices that are on a chromsome arm, and are thus
	 * block diagonal. This is assumed to be true and is not checked.
	 *
	 * @param genes map as created by IoUtils.readGenesAsChrArmMap
	 * @param rowOrder A list with gene ids in the order of the data.
	 * @return
	 * @throws IOException
	 */
	public static List<int[]> createBlockDiagonalIndexFromGenes(Map<String, List<Gene>> genes, List<String> rowOrder) throws IOException {

		// Determine the indices in the matrix that lie together on a chromosome arm
		List<int[]> blockDiagonalIndices = new ArrayList<>();
		Map<String, TIntArrayList> blockDiagonalIndicesTmp = new HashMap<>();

        LOGGER.info("Read genes on " + genes.size() + " blocks");

        int totalCount = 0;

		for (String key : genes.keySet()) {
			for (Gene curGene : genes.get(key)) {
				int idx = rowOrder.indexOf(curGene.getGene());
				if (idx > -1) {
					if (!blockDiagonalIndicesTmp.containsKey(key)) {
						blockDiagonalIndicesTmp.put(key, new TIntArrayList());
					}
					blockDiagonalIndicesTmp.get(key).add(idx);
					totalCount++;
				}
			}
		}

		if (rowOrder.size() > totalCount) {
			throw new IllegalArgumentException("There were genes on rows not found in --genes file. "
					+ "Please make sure all genes are available in --genes or use -ro to make select a subset.");
		}

		// Convert to int[] so it's easier to use with DoubleMatrix2D;
		for (TIntArrayList curList : blockDiagonalIndicesTmp.values()) {
			blockDiagonalIndices.add(curList.toArray());
		}

        LOGGER.info("Created index for " + totalCount + " genes on " + blockDiagonalIndices.size() + " blocks");

        return blockDiagonalIndices;
    }


    /**
     * Perform eigen decomposition of a block diagonal matrix. Will decompose each block seperately according to the
     * blocks defined in index. Can use either Pcolt or Jblas for eigen decomposition.
     * @param provider
     * @param index
     * @param useJblas
     * @return
     * @throws UnsatisfiedLinkError
     */
    public static DoubleMatrixDataset<String, String>[] blockDiagonalEigenDecomposition(BlockDiagonalDoubleMatrixProvider provider, List<int[]> index, boolean useJblas) throws UnsatisfiedLinkError, Exception {

        List<IndexedDouble> eigenvalues = new ArrayList<>(provider.columns());
        List<String> eigenvectorNames = new ArrayList<>(provider.columns());
        for (int i = 0; i < provider.columns(); i++) {
            eigenvectorNames.add("V" + i);
        }
        List<String> eigenvalueNames = new ArrayList<>(1);
        eigenvalueNames.add("eigenvalues");

        DoubleMatrixDataset<String, String> U = new DoubleMatrixDataset<>(provider.getRowNames(), eigenvectorNames);
        DoubleMatrixDataset<String, String> L = new DoubleMatrixDataset<>(eigenvectorNames, eigenvalueNames);

        ProgressBar pb = new ProgressBar("Eigen decomposition per block", index.size());

        int masterIndex = 0;

        if (useJblas) {
            // Use Jblas for eigen decompotision
            try {
                for (int curBlock = 0; curBlock < index.size(); curBlock++) {
                    int[] curIndex = index.get(curBlock);
                    //DoubleMatrix curMatrix = new DoubleMatrix(sigma.getMatrix().viewSelection(curIndex, curIndex).toArray());
                    DoubleMatrix curMatrix = toJblasDoubleMatrix(provider.viewBlock(curIndex));
                    DoubleMatrix[] eigen = Eigen.symmetricEigenvectors(curMatrix);

                    for (int i = 0; i < curMatrix.rows; i++) {
                        eigenvalues.add(new IndexedDouble(eigen[1].get(i, i), masterIndex));
                        masterIndex++;
                        for (int j = 0; j < curMatrix.columns; j++) {
                            U.setElementQuick(curIndex[i], curIndex[j], eigen[0].get(i, j));
                        }
                    }
                    pb.step();
                }

            } catch (UnsatisfiedLinkError e) {
                pb.close();
                throw e;
            }

        } else {
            // Use P-colt for eigen decomposition
            for (int curBlock = 0; curBlock < index.size(); curBlock++) {
                int[] curIndex = index.get(curBlock);
                DoubleMatrix2D curMatrix = provider.viewBlock(curIndex);
                DenseDoubleEigenvalueDecomposition eigen = new DenseDoubleEigenvalueDecomposition(curMatrix);

                for (int i = 0; i < curMatrix.rows(); i++) {
                    eigenvalues.add(new IndexedDouble(eigen.getRealEigenvalues().get(i), masterIndex));
                    masterIndex++;
                    for (int j = 0; j < curMatrix.columns(); j++) {
                        U.setElementQuick(curIndex[i], curIndex[j], eigen.getV().getQuick(i, j));
                    }
                }

                pb.step();
            }
        }

        // Order according to eigenvalues, large to small
        orderToEigenvalues(eigenvalues, U, L);

        pb.close();

        DoubleMatrixDataset<String, String>[] output = new DoubleMatrixDataset[2];
        output[0] = L;
        output[1] = U;
        return output;
    }


    /**
     * Eigen decomposition of DoubleMatrixDataset. Ensures the order of eigenvalues and vectors is decreasing in value.
     * @param Sigma
     * @param useJblas
     * @return
     * @throws Exception
     */
    public static DoubleMatrixDataset<String, String>[] eigenDecomposition(DoubleMatrixDataset<String, String> Sigma, boolean useJblas) throws Exception {
        LOGGER.warn("Performing eigen decomposition on full sigma matrix. Is this what you want?");

        List<IndexedDouble> eigenvalues = new ArrayList<>(Sigma.columns());
        List<String> eigenvectorNames = new ArrayList<>(Sigma.columns());
        for (int i = 0; i < Sigma.columns(); i++) {
            eigenvectorNames.add("V" + i);
        }
        List<String> eigenvalueNames = new ArrayList<>(1);
        eigenvalueNames.add("eigenvalues");

        DoubleMatrixDataset<String, String> U = new DoubleMatrixDataset<>(Sigma.getRowObjects(), eigenvectorNames);
        DoubleMatrixDataset<String, String> L = new DoubleMatrixDataset<>(eigenvectorNames, eigenvalueNames);

        if (useJblas) {
            DoubleMatrix curMatrix = toJblasDoubleMatrix(Sigma);
            DoubleMatrix[] eigen = Eigen.symmetricEigenvectors(curMatrix);

            for (int i = 0; i < curMatrix.rows; i++) {
                eigenvalues.add(new IndexedDouble(eigen[1].get(i, i), i));
                for (int j = 0; j < curMatrix.columns; j++) {
                    U.setElementQuick(i, j, eigen[0].get(i, j));
                }
            }

        } else {
            DenseDoubleEigenvalueDecomposition eigen = new DenseDoubleEigenvalueDecomposition(Sigma.getMatrix());
            U.setMatrix(eigen.getV());
            for (int i=0; i<eigen.getRealEigenvalues().size(); i++) {
                eigenvalues.add(new IndexedDouble(eigen.getRealEigenvalues().get(i), i));
            }
        }

        // Order according to eigenvalues, large to small
        orderToEigenvalues(eigenvalues, U, L);

        DoubleMatrixDataset<String, String>[] output = new DoubleMatrixDataset[2];
        output[0] = L;
        output[1] = U;

        return output;
    }


    /**
     * Orders a eigen decomposition to the eigenvalues in decreasing order.
     * @param eigenvalues   List of eigenvalues as indexed doubles
     * @param U
     * @param L
     * @throws Exception
     */
    private static void orderToEigenvalues(List<IndexedDouble> eigenvalues, DoubleMatrixDataset<String, String> U, DoubleMatrixDataset<String, String> L) throws Exception {

        // Keep track of original eigenvector names
        List<String> eigenvectorNames = new ArrayList<>(U.getColObjects());

        // Order according to eigenvalues, large to small
        eigenvalues.sort(Comparator.comparingDouble(IndexedDouble::getValue));
        Collections.reverse(eigenvalues);

        LinkedHashMap<String, Integer> eigenvectorOrder = new LinkedHashMap<>();
        int i = 0;
        for (IndexedDouble cur : eigenvalues) {
            eigenvectorOrder.put(eigenvectorNames.get(cur.getIndex()), i);
            L.setElementQuick(i, 0, cur.getValue());
            i++;
        }

        U.reorderCols(eigenvectorOrder);

        // Now that eigenvectors are properly ordered, reset the names.
        U.setColObjects(eigenvectorNames);
    }


    /**
     * Convert a DoubleMatrixDataset to Jblas DoubleMatrix
     * @param A
     * @return
     */
    private static DoubleMatrix toJblasDoubleMatrix(DoubleMatrixDataset A) {
        DoubleMatrix output = new DoubleMatrix(A.rows(), A.columns());

        for (int i=0; i< A.rows(); i++) {
            for (int j = 0; j < A.columns(); j++) {
                output.put(i, j, A.getElementQuick(i, j));
            }
        }
        return output;
    }

    /**
     * Convert a colt DoubleMatrix2d to Jblas DoubleMatrix
     * @param A
     * @return
     */
    private static DoubleMatrix toJblasDoubleMatrix(DoubleMatrix2D A) {
        DoubleMatrix output = new DoubleMatrix(A.rows(), A.columns());

        for (int i=0; i< A.rows(); i++) {
            for (int j = 0; j < A.columns(); j++) {
                output.put(i, j, A.getQuick(i, j));
            }
        }
        return output;
    }


    /**
     * Calculate A %*% B for a matrix A whose structure has some form of block diagonality. I.e. there are rows
     * and columns that are exactly zero. @param index gives the indices to multiply together.
     * JBlas is much quicker on larger matrices but has memory overhead compared to colt matrix mult.
     * <p>
     * Columns in A are first subset per block. The rows to keep in A are determined by removing rows that have only
     * zero values in the  subset of A's columns.
     * <p>
     * Indices provided in index should match the columns of A and the rows of B.
     *
     * @param A
     * @param B
     * @param index
     * @param useJblas
     * @return
     */
    private static DoubleMatrix2D blockDiagonalMult(DoubleMatrix2D A, DoubleMatrix2D B, List<int[]> index, boolean useJblas) {
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
     * <p>
     * Columns in A are first subset per block. The rows to keep in A are determined by removing rows that have only
     * zero values in the  subset of A's columns.
     * <p>
     * Indices provided in index should match the columns of A and the rows of B.
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
            //Acache[i] = new DoubleMatrix(Ac.viewSelection(null, index.get(i)).toArray());
            Acache[i] = toJblasDoubleMatrix(Ac.viewSelection(null, index.get(i)));
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
     * <p>
     * Columns in A are first subset per block. The rows to keep in A are determined by removing rows that have only
     * zero values in the  subset of A's columns.
     * <p>
     * Indices provided in index should match the columns of A and the rows of B.
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
	 * Calculate the pairwise product between a matrix and a "diagonal" matrix
	 * column wise. The diagonal is represented as a vector.
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
	 * Calculate the pairwise product between a vector and a diagonal matrix
	 * column wise. The diagonal is represented as a vector.
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
