package nl.systemsgenetics.downstreamer.runners;

import cern.colt.GenericPermuting;
import cern.colt.function.tdouble.DoubleDoubleFunction;
import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import cern.colt.matrix.tdouble.algo.decomposition.DenseDoubleEigenvalueDecomposition;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import cern.jet.math.tdouble.DoubleFunctions;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TObjectIntHashMap;
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
import java.util.logging.Level;
import java.util.stream.IntStream;
import me.tongfei.progressbar.ProgressBarStyle;
import static nl.systemsgenetics.downstreamer.runners.DownstreamerEnrichment.inplacePvalueToZscore;
import nl.systemsgenetics.downstreamer.runners.options.OptionsBase;

public class DownstreamerRegressionEngine {

	private static final Logger LOGGER = LogManager.getLogger(DownstreamerRegressionEngine.class);

	public static final String MAIN_EFFECT_COL_NAME = "main_effect";
	public static final String INTERCEPT_COL_NAME = "intercept";

	static File debugFolder;

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

		inplacePvalueToZscore(Y);

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
		LinkedHashMap<String, ArrayList<String>> blockDiagonalIndices2 = null;

		if (options.hasGenes()) {
			blockDiagonalIndices = createBlockDiagonalIndexFromGenes(options.getGenes(), finalRowSelection);
			blockDiagonalIndices2 = createBlockDiagonalIndexFromGenes2(options.getGenes(), finalRowSelection);

//			if (options.isDebugMode()) {
//				IoUtils.writeBlockDiagonalIndices(blockDiagonalIndices, finalRowSelection, options.getDebugFolder() + "/block_diagonal_indices.txt");
//			}
		}

		// Depending on input, first decompse Sigma, or run with precomputed eigen decomp.
		List<LinearRegressionResult> output;

		// TODO: in future make sure this doesn't read the full matrix
		BlockDiagonalDoubleMatrixProvider sigmaProvider = new DoubleMatrixDatasetBlockDiagonalProvider(Sigma);
		DoubleMatrixDataset<String, String>[] eigen = null;

		if (Sigma != null && blockDiagonalIndices != null) {
			eigen = blockDiagonalEigenDecomposition(finalRowSelection, sigmaProvider, blockDiagonalIndices2, options.useJblas());
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
				for (String cur : finalRowSelection) {
					writer.write(cur);
					writer.newLine();
				}
				writer.flush();
				writer.close();
			}
		}

		output = performDownstreamerRegression(X, Y, C, U, L,
				blockDiagonalIndices2,
				options.getPercentageOfVariance(),
				options.fitIntercept(),
				options.regressCovariates(),
				options.useJblas(), 0);

		// Save linear regression results
		for (LinearRegressionResult curRes : output) {
			curRes.save(options.getOutputBasePath() + "_" + curRes.getName(), false);
		}
	}

	/**
	 * Runs DS regression if eigendecompostion has been pre-computed.These are
	 * given by U and L.Expects that these still need to be filtered. If they
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
	 * @param permutations use 0 for no permutations and >0 to permute Y in
	 * eigenvector space.
	 * @return Per trait (column in Y) a LinearRegressionResult.
	 */
	public static List<LinearRegressionResult> performDownstreamerRegression(DoubleMatrixDataset<String, String> X,
			DoubleMatrixDataset<String, String> Y,
			DoubleMatrixDataset<String, String> C,
			DoubleMatrixDataset<String, String> U,
			DoubleMatrixDataset<String, String> L,
			final LinkedHashMap<String, ArrayList<String>> blockDiagonalIndices,
			double percentageOfVariance,
			boolean fitIntercept,
			boolean regressCovariates,
			boolean useJblas,
			int permutations) {

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

		return performDownstreamerRegression(X, Y, C, U, L, blockDiagonalIndices, fitIntercept, regressCovariates, useJblas, permutations);
	}

	/**
	 * Runs DS regression if U and L have been pre-filtered on only the
	 * eigenvectors to use. See above for parameter descriptions as they are
	 * identical.
	 */
	public static List<LinearRegressionResult> performDownstreamerRegression(final DoubleMatrixDataset<String, String> X,
			final DoubleMatrixDataset<String, String> Y,
			final DoubleMatrixDataset<String, String> C,
			final DoubleMatrixDataset<String, String> U,
			final DoubleMatrixDataset<String, String> L,
			final LinkedHashMap<String, ArrayList<String>> blockDiagonalIndices,
			final boolean fitIntercept,
			final boolean regressCovariates,
			final boolean useJblas,
			int permutations) {

		debugFolder = OptionsBase.getDebugFolder();
//
//		try {
//			X.save(new File(debugFolder, "X.txt"));
//			Y.save(new File(debugFolder, "Y.txt"));
//			if (C != null) {
//				C.save(new File(debugFolder, "C.txt"));
//			}
//			U.save(new File(debugFolder, "U.txt"));
//			L.save(new File(debugFolder, "L.txt"));
//		} catch (IOException ex) {
//			throw new RuntimeException();
//		}

		if (X.rows() != Y.rows()) {
			throw new RuntimeException("Internal error: X and Y should contain equal number of rows");
		}
		if (C != null && C.rows() != Y.rows()) {
			throw new RuntimeException("Internal error: C and Y should contain equal number of rows");
		}
		if (U.rows() != Y.rows()) {
			throw new RuntimeException("Internal error: U and Y should contain equal number of rows");
		}
		if (U.columns() != L.rows()) {
			throw new RuntimeException("Internal error: U columns and L rows should be equal");
		}

		if (permutations > 0 && Y.columns() > 1) {
			throw new RuntimeException("Max 1 trait per permutations round (this should never happen in normal usage)");
		}

//		if (!X.getHashRows().keySet().equals(Y.getHashRows().keySet())) {
//			throw new RuntimeException("Internal error: X and Y should contain row elements in same order");
//		}
//		if (!X.getHashRows().keySet().equals(Y.getHashRows().keySet())) {
//			throw new RuntimeException("Internal error: C and Y should contain row elements in same order");
//		}
//		if (!U.getHashRows().keySet().equals(Y.getHashRows().keySet())) {
//			throw new RuntimeException("Internal error: U and Y should contain row elements in same order");
//		}
//		if (!U.getHashCols().keySet().equals(L.getHashRows().keySet())) {
//			throw new RuntimeException("Internal error: U columns and L rows should contain elements in same order");
//		}
		// Pre-transpose U as this can be re-used as well
		DoubleMatrixDataset<String, String> UHatT = U.viewDice();
		//UHatT = DoubleFactory2D.sparse.make(UHatT.toArray());

		// Calculate inverse of eigenvalues
		DoubleMatrix1D LHatInv = L.viewCol(0).copy();
		LHatInv.assign(DoubleFunctions.inv);

		logInfoMem("Starting regression for " + X.columns() + " pathways.");

		// Determine the degrees of freedom. -1 for main_effect
		int degreesOfFreedomTmp = UHatT.rows() - 1;

		// Set the names of the predictors in the correct order
		// and determine the degrees of freedom
		List<String> predictorNames = new ArrayList<>();
		if (fitIntercept) {
			predictorNames.add(INTERCEPT_COL_NAME);
			degreesOfFreedomTmp = degreesOfFreedomTmp - 1;
		}
		predictorNames.add(MAIN_EFFECT_COL_NAME);
		if (C != null) {
			predictorNames.addAll(C.getColObjects());
			degreesOfFreedomTmp = degreesOfFreedomTmp - C.columns();
		}

		final int degreesOfFreedom = degreesOfFreedomTmp;

		// Pre calculate predictors and covariates in eigenvectors space.
		final DoubleMatrixDataset XhatCache;
		final DoubleMatrixDataset ChatCache;

		// Pre-calculate XHat and Xhat in a block diagonal way, or using full matrix mult.
		if (blockDiagonalIndices != null) {
			logInfoMem("Precomputing XHat and CHat");

			XhatCache = blockDiagonalMult(UHatT, X, blockDiagonalIndices, useJblas);

			if (C != null) {
				ChatCache = blockDiagonalMult(UHatT, C, blockDiagonalIndices, useJblas);
			} else {
				ChatCache = null;
			}
			logInfoMem("Done");
		} else {
			XhatCache = null;
			if (C != null) {
				ChatCache = new DoubleMatrixDataset(mult(UHatT.getMatrix(), C.getMatrix()));
			} else {
				ChatCache = null;
			}
		}

//		try {
//			if (XhatCache != null) {
//				XhatCache.save(new File(debugFolder, "XhatCache.txt"));
//			}
//		} catch (Exception ex) {
//			throw new RuntimeException(ex);
//		}
		final int nrTraits;
		final int[] permutationInts;
		if (permutations == 0) {
			nrTraits = Y.columns();
			permutationInts = null;
		} else {
			Random r = new Random(42);
			permutationInts = new int[permutations];
			for (int p = 0; p < permutations; ++p) {
				permutationInts[p] = Math.abs(r.nextInt()) + 1;
			}
			nrTraits = permutations;
		}

		ProgressBar pb = new ProgressBar("Linear regressions", X.columns() * nrTraits, ProgressBarStyle.ASCII);

		// Store output
		final LinearRegressionResult[] resultsArray = new LinearRegressionResult[nrTraits];

		IntStream.range(0, nrTraits).parallel().forEach(curY -> {
			//for (int curY = 0; curY < Y.columns(); curY++) {

			// Keep track of current covariates, if these are used to regress need to keep ChatCache.
			// but curCovariates can be set to null
			DoubleMatrixDataset<String, String> curCovariates = ChatCache;
			List<String> curPredictorNames = predictorNames;

			// Pre-compute Y^ as this can be re-used
			final DoubleMatrix2D y;
			final String name;
			if (permutations == 0) {
				y = Y.viewColAsMmatrix(curY);
				name = Y.getColObjects().get(curY);
			} else {
				y = Y.viewColAsMmatrix(0);
				y.viewColumn(0).assign(y.viewColumn(0).viewSelection(GenericPermuting.permutation(permutationInts[curY], Y.rows())));
				name = "P" + nrTraits;
			}
			final DoubleMatrix2D YHat = mult(UHatT.getMatrix(), y);

//			if (curY == 0) {
//				try {
//					new DoubleMatrixDataset<String, String>(YHat).save(new File(debugFolder, "YHat.txt"));
//				} catch (IOException ex) {
//					throw new RuntimeException(ex);
//				}
//			}
			// Instead of including covariates in the model, first regress their effects
			if (C != null && regressCovariates) {

				inplaceDownstreamerRegressionResidualsPrecomp(YHat, ChatCache.getMatrix(), LHatInv);

				// Set covariates to null so they are not used for the rest of the regression
				curCovariates = null;
				curPredictorNames = predictorNames.subList(0, predictorNames.size() - C.columns());
			}

			// Object to save output
			final LinearRegressionResult result = new LinearRegressionResult(X.getColObjects(),
					curPredictorNames,
					degreesOfFreedom,
					name,
					UHatT.rows()
			);

			// Calculate beta's and SE for each pathway
			for (int curPathway = 0; curPathway < X.columns(); curPathway++) {

				// Get the current pathway, if not pre-multiplied by eigenvectors, do so now.
				DoubleMatrix2D XCur;
				if (XhatCache == null) {
					XCur = mult(UHatT.getMatrix(), X.viewColAsMmatrix(curPathway));
				} else {
					XCur = XhatCache.getMatrix().viewSelection(null, new int[]{curPathway});
				}

				// If covariates are provided add them to the design matrix
				if (curCovariates != null) {
					XCur = DoubleFactory2D.dense.appendColumns(XCur, ChatCache.getMatrix());
				}

//				if (curPathway == 0) {
//					try {
//						new DoubleMatrixDataset(XCur).save(new File(debugFolder, "XCur.txt"));
//					} catch (IOException ex) {
//						throw new RuntimeException(ex);
//					}
//				}
				double[] curRes;
				try {

					curRes = downstreamerRegressionPrecomp(XCur,
							YHat,
							LHatInv,
							fitIntercept);

				} catch (IOException ex) {
					throw new RuntimeException(ex);
				}

				result.appendBetas(curPathway, ArrayUtils.subarray(curRes, 0, curRes.length / 2));
				result.appendSe(curPathway, ArrayUtils.subarray(curRes, curRes.length / 2, curRes.length));

				pb.step();
			}

//			try {
//
//				DoubleMatrixDataset<String, String> betas = result.getBeta();
//
//				betas.save(new File(debugFolder, "res_beta.txt"));
//				result.getStandardError().save(new File(debugFolder, "res_se.txt"));
//
//				DoubleMatrixDataset t = new DoubleMatrixDataset(betas.getRowObjects(), Arrays.asList(new String[]{"T"}));
//				t.getMatrix().viewColumn(0).assign(result.getTstatForMainEffect());
//				t.save(new File(debugFolder, "res_t.txt"));
//
//			} catch (IOException ex) {
//				throw new RuntimeException();
//			}

			resultsArray[curY] = result;
			//results.add(result);
		});

		pb.close();
		logInfoMem("Done with regression for " + nrTraits + " traits and " + X.columns() + " pathways.");

		return Arrays.asList(resultsArray);
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
	private static double[] downstreamerRegressionPrecomp(DoubleMatrix2D XHat, DoubleMatrix2D YHat, DoubleMatrix1D LHatInv, boolean fitIntercept) throws IOException {

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
//		new DoubleMatrixDataset(YHat).save(new File(debugFolder, "YHat_2.txt"));
//		new DoubleMatrixDataset(LHatInv.reshape(1, (int) LHatInv.size())).save(new File(debugFolder, "LHatInv_2.txt"));
//		new DoubleMatrixDataset(XHat).save(new File(debugFolder, "XHat_2.txt"));
		DoubleMatrix2D a = mult(transpose(multDiag(YHat, LHatInv)), XHat);
//		

//		new DoubleMatrixDataset(a).save(new File(debugFolder, "a.txt"));
//		LOGGER.debug("XHat " + XHat.toStringShort());
//		LOGGER.debug(XHat.toString());
//		LOGGER.debug("LHatInv " + LHatInv.toStringShort());
//		
//		
//		LOGGER.debug("Test");
//		DoubleMatrix2D test = mult(transpose(multDiag(XHat, LHatInv)), XHat);
//		LOGGER.debug(test.toString());
//		
		DoubleMatrix2D b = inverse(mult(transpose(multDiag(XHat, LHatInv)), XHat));

//		new DoubleMatrixDataset(b).save(new File(debugFolder, "b.txt"));
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

//		LOGGER.debug("Y " + Y.toString());
//		LOGGER.debug("X " + X.toString());
//		LOGGER.debug("LHatInv " + LHatInv.toString());
		//DoubleMatrix2D design = DoubleFactory2D.dense.appendColumns(DoubleFactory2D.dense.make(X.rows(), 1, 1), X);
		DoubleMatrix2D design = X.copy();

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
		Map<String, TIntArrayList> blockDiagonalIndicesTmp = new LinkedHashMap<>();

		LOGGER.debug("Read genes on " + genes.size() + " blocks");

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

		LOGGER.debug("Created index for " + totalCount + " genes on " + blockDiagonalIndices.size() + " blocks");

		return blockDiagonalIndices;
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
	public static LinkedHashMap<String, ArrayList<String>> createBlockDiagonalIndexFromGenes2(File geneFile, List<String> rowOrder) throws IOException {
		Map<String, List<Gene>> genes = IoUtils.readGenesAsChrArmMap(geneFile);
		return createBlockDiagonalIndexFromGenes2(genes, rowOrder);
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
	public static LinkedHashMap<String, ArrayList<String>> createBlockDiagonalIndexFromGenes2(Map<String, List<Gene>> genes, List<String> rowOrder) throws IOException {

		LinkedHashMap<String, ArrayList<String>> blockDiagonalIndices2 = new LinkedHashMap<>();

		LOGGER.debug("Read genes on " + genes.size() + " blocks");

		int totalCount = 0;

		for (String key : genes.keySet()) {

			HashSet<String> armGeneSet = new HashSet<>();
			for (Gene curGene : genes.get(key)) {
				armGeneSet.add(curGene.getGene());
			}

			for (String geneName : rowOrder) {
				if (armGeneSet.contains(geneName)) {

					if (!blockDiagonalIndices2.containsKey(key)) {
						blockDiagonalIndices2.put(key, new ArrayList<String>());
					}
					blockDiagonalIndices2.get(key).add(geneName);
					totalCount++;

				}
			}

		}

		if (rowOrder.size() > totalCount) {
			throw new IllegalArgumentException("There were genes on rows not found in --genes file. "
					+ "Please make sure all genes are available in --genes or use -ro to make select a subset.");
		}

		LOGGER.debug("Created index for " + totalCount + " genes on " + blockDiagonalIndices2.size() + " blocks");

		return blockDiagonalIndices2;
	}

	/**
	 * Perform eigen decomposition of a block diagonal matrix.Will decompose
	 * each block seperately according to the blocks defined in index. Can use
	 * either Pcolt or Jblas for eigen decomposition.
	 *
	 * @param allGenes
	 * @param provider
	 * @param index
	 * @param useJblas
	 * @return [0] eigenvalues [1] eigenvectors
	 * @throws UnsatisfiedLinkError
	 */
	public static DoubleMatrixDataset<String, String>[] blockDiagonalEigenDecomposition(List<String> allGenes, BlockDiagonalDoubleMatrixProvider provider, LinkedHashMap<String, ArrayList<String>> index, boolean useJblas) throws UnsatisfiedLinkError, Exception {

		final IndexedDouble[] eigenvalues = new IndexedDouble[allGenes.size()];
		List<String> eigenvectorNames = new ArrayList<>(allGenes.size());
		for (int i = 1; i <= allGenes.size(); i++) {
			eigenvectorNames.add("V" + i);
		}
		List<String> eigenvalueNames = new ArrayList<>(1);
		eigenvalueNames.add("eigenvalues");

		final DoubleMatrixDataset<String, String> U = new DoubleMatrixDataset<>(allGenes, eigenvectorNames);
		final DoubleMatrixDataset<String, String> L = new DoubleMatrixDataset<>(eigenvectorNames, eigenvalueNames);

		ProgressBar pb = new ProgressBar("Eigen decomposition per chromosome arm", index.size(), ProgressBarStyle.ASCII);

		TObjectIntMap<String> columnIndexStartOfBlock = new TObjectIntHashMap<>(index.size());
		int masterIndex = 0;
		for (Map.Entry<String, ArrayList<String>> block : index.entrySet()) {
			columnIndexStartOfBlock.put(block.getKey(), masterIndex);
			//LOGGER.debug("block" + block.getKey() + " genes: " + block.getValue().size() + " index: " + masterIndex);
			masterIndex += block.getValue().size();

		}

		if (useJblas) {
			// Use Jblas for eigen decompotision
			try {
				//NOTE: this can't be done parralel because Eigen.symmetricEigenvectors uses shared static variables
				for (Map.Entry<String, ArrayList<String>> block : index.entrySet()) {

					int colIndex = columnIndexStartOfBlock.get(block.getKey());

					final DoubleMatrix curMatrix = toJblasDoubleMatrix(provider.viewBlock(block.getKey(), block.getValue()).getMatrix());
					final DoubleMatrix[] eigen = Eigen.symmetricEigenvectors(curMatrix);

					//get U in the same order 
					final DoubleMatrixDataset<String, String> UBlock = U.viewRowSelection(block.getValue());

					for (int j = 0; j < curMatrix.columns; j++) {
						eigenvalues[colIndex] = new IndexedDouble(eigen[1].get(j, j), colIndex);
						for (int i = 0; i < curMatrix.rows; i++) {
							UBlock.setElementQuick(i, colIndex, eigen[0].get(i, j));
						}
						colIndex++;
					}

					pb.step();
				}

			} catch (UnsatisfiedLinkError e) {
				pb.close();
				throw e;
			}

		} else {
			// Use P-colt for eigen decomposition

			for (Map.Entry<String, ArrayList<String>> block : index.entrySet()) {
			//index.entrySet().parallelStream().forEach((Map.Entry<String, ArrayList<String>> block) -> {
				try {

					int colIndex = columnIndexStartOfBlock.get(block.getKey());

					//for (Map.Entry<String, ArrayList<String>> block : index.entrySet()) {
					//for (int curBlock = 0; curBlock < index.size(); curBlock++) {
					final DoubleMatrix2D curMatrix = provider.viewBlock(block.getKey(), block.getValue()).getMatrix().copy();

					//this functions seems to not contain multithreading
					final DenseDoubleEigenvalueDecomposition eigen;
					try {
						eigen = new DenseDoubleEigenvalueDecomposition(curMatrix);
					} catch (RuntimeException ex) {
						System.out.println(OptionsBase.getDebugFolder()  + "test.txt");
						provider.viewBlock(block.getKey(), block.getValue()).save(OptionsBase.getDebugFolder()  + "test.txt");
						System.out.println("curMatrix: " + curMatrix.toStringShort());
						System.out.println(curMatrix.getQuick(0, 0));
						System.out.println(curMatrix.getQuick(0, 1));
						System.out.println(curMatrix.getQuick(curMatrix.rows()-1, curMatrix.columns()-1));
						System.out.println(curMatrix.getQuick(curMatrix.rows()-1, curMatrix.columns()-2));
						System.out.println("test");
						throw new Exception(ex);
					}

					//get U in the same order 
					DoubleMatrixDataset<String, String> UBlock = U.viewRowSelection(block.getValue());

					for (int j = 0; j < curMatrix.columns(); j++) {
						eigenvalues[colIndex] = new IndexedDouble(eigen.getRealEigenvalues().get(j), colIndex);
						for (int i = 0; i < curMatrix.rows(); i++) {
							UBlock.setElementQuick(i, colIndex, eigen.getV().getQuick(i, j));
						}
						colIndex++;
					}

					pb.step();
				} catch (Exception ex) {
					pb.close();
					throw (new RuntimeException(ex));
				}
			}//);
		}

		//U.printMatrix();
		// Order according to eigenvalues, large to small. eignevalues and L are inplace U is returned as a view
		DoubleMatrixDataset<String, String> U2 = orderToEigenvalues(Arrays.asList(eigenvalues), U, L);

		pb.close();

		DoubleMatrixDataset<String, String>[] output = new DoubleMatrixDataset[2];
		output[0] = L;
		output[1] = U2;
		return output;
	}

	/**
	 * Eigen decomposition of DoubleMatrixDataset. Ensures the order of
	 * eigenvalues and vectors is decreasing in value.
	 *
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
			for (int i = 0; i < eigen.getRealEigenvalues().size(); i++) {
				eigenvalues.add(new IndexedDouble(eigen.getRealEigenvalues().get(i), i));
			}
		}

		// Order according to eigenvalues, large to small
		U = orderToEigenvalues(eigenvalues, U, L);

		DoubleMatrixDataset<String, String>[] output = new DoubleMatrixDataset[2];
		output[0] = L;
		output[1] = U;

		return output;
	}

	/**
	 * Orders a eigen decomposition to the eigenvalues in decreasing order.
	 *
	 * @param eigenvalues List of eigenvalues as indexed doubles
	 * @param U
	 * @param L
	 * @throws Exception
	 */
	private static DoubleMatrixDataset<String, String> orderToEigenvalues(List<IndexedDouble> eigenvalues, DoubleMatrixDataset<String, String> U, DoubleMatrixDataset<String, String> L) throws Exception {

//		System.out.println("Ü");
//		U.printSummary();
//		
//		System.out.println("L");
//		L.printSummary();
//		
		// Keep track of original eigenvector names
		List<String> eigenvectorNames = new ArrayList<>(U.getColObjects());

		// Order according to eigenvalues, large to small
		eigenvalues.sort(Comparator.comparingDouble(IndexedDouble::getValue));
		Collections.reverse(eigenvalues);

		ArrayList<String> eigenvectorOrder = new ArrayList<>(eigenvalues.size());
		int i = 0;
		for (IndexedDouble cur : eigenvalues) {
			eigenvectorOrder.add(eigenvectorNames.get(cur.getIndex()));
			L.setElementQuick(i, 0, cur.getValue());
			i++;
		}

		DoubleMatrixDataset<String, String> U2 = U.viewColSelection(eigenvectorOrder);
		U2.setColObjects(eigenvectorNames);

//		System.out.println("Ü2");
//		U2.printSummary();
//		
//		
//		System.out.println("L2");
//		L.printSummary();
//		
		return U2;

	}

	/**
	 * Convert a DoubleMatrixDataset to Jblas DoubleMatrix
	 *
	 * @param A
	 * @return
	 */
	private static DoubleMatrix toJblasDoubleMatrix(DoubleMatrixDataset A) {
		return (toJblasDoubleMatrix(A.getMatrix()));
	}

	/**
	 * Convert a colt DoubleMatrix2d to Jblas DoubleMatrix
	 *
	 * @param A
	 * @return
	 */
	private static DoubleMatrix toJblasDoubleMatrix(DoubleMatrix2D A) {
		DoubleMatrix output = new DoubleMatrix(A.rows(), A.columns());

		for (int i = 0; i < A.rows(); i++) {
			for (int j = 0; j < A.columns(); j++) {
				output.put(i, j, A.getQuick(i, j));
			}
		}
		return output;
	}

	/**
	 * Calculate A %*% B for a matrix A whose structure has some form of block
	 * diagonality. I.e. there are rows and columns that are exactly zero.
	 *
	 * @param index gives the indices to multiply together. JBlas is much
	 * quicker on larger matrices but has memory overhead compared to colt
	 * matrix mult.
	 * <p>
	 * Columns in A are first subset per block. The rows to keep in A are
	 * determined by removing rows that have only zero values in the subset of
	 * A's columns.
	 * <p>
	 * Indices provided in index should match the columns of A and the rows of
	 * B.
	 *
	 * @param A
	 * @param B
	 * @param index
	 * @param useJblas
	 * @return
	 */
	private static DoubleMatrixDataset<String, String> blockDiagonalMult(DoubleMatrixDataset<String, String> A, DoubleMatrixDataset<String, String> B, final LinkedHashMap<String, ArrayList<String>> blockDiagonalIndices, boolean useJblas) {
		if (useJblas) {
			return multBlockDiagonalSubsetPerColJblas(A, B, blockDiagonalIndices);
		} else {
			return multBlockDiagonalSubsetPerColPcolt(A, B, blockDiagonalIndices);
		}
	}

	/**
	 * Calculate A %*% B for a matrix A whose structure has some form of block
	 * diagonality. I.e. there are rows and columns that are exactly zero.
	 *
	 * @param index gives the indices to multiply together. Uses Jblas for
	 * matrix multiplications. This does give memory overhead but is much
	 * quciker.
	 * <p>
	 * Columns in A are first subset per block. The rows to keep in A are
	 * determined by removing rows that have only zero values in the subset of
	 * A's columns.
	 * <p>
	 * Indices provided in index should match the columns of A and the rows of
	 * B.
	 *
	 * @param Ac
	 * @param Bc
	 * @param index
	 * @return
	 */
	private static DoubleMatrixDataset<String, String> multBlockDiagonalSubsetPerColJblas(final DoubleMatrixDataset<String, String> Ac, final DoubleMatrixDataset<String, String> Bc, final LinkedHashMap<String, ArrayList<String>> blockDiagonalIndices) {

		final DoubleMatrixDataset C = new DoubleMatrixDataset(Ac.getRowObjects(), Bc.getColObjects());//DoubleFactory2D.dense.make(Ac.rows(), Bc.columns());

		final int nrBlocks = blockDiagonalIndices.size();
		ProgressBar pb = new ProgressBar("Multiplying per block", Bc.columns() * nrBlocks, ProgressBarStyle.ASCII);

		final DoubleMatrix[] Acache = new DoubleMatrix[nrBlocks];
		final DoubleMatrix[] Bcache = new DoubleMatrix[nrBlocks];

		int[][] AcacheIndices = new int[nrBlocks][];

		final ArrayList<String> blocks = new ArrayList(blockDiagonalIndices.keySet());

		// Compile the subset caches. Precomputing saves overhead, but is more memory intense.
		for (int i = 0; i < nrBlocks; i++) {

			final ArrayList<String> blockGenes = blockDiagonalIndices.get(blocks.get(i));

			//Acache[i] = new DoubleMatrix(Ac.viewSelection(null, index.get(i)).toArray());
			Acache[i] = toJblasDoubleMatrix(Ac.viewColSelection(blockGenes));
			DoubleMatrix tmp = Acache[i].rowSums();
			ArrayList<Integer> tmpIndices = new ArrayList<>();
			for (int j = 0; j < tmp.rows; j++) {
				if (tmp.get(j, 0) != 0) {
					tmpIndices.add(j);
				}
			}

			AcacheIndices[i] = tmpIndices.stream().mapToInt(p -> p).toArray();
			Acache[i] = Acache[i].getRows(AcacheIndices[i]);
			Bcache[i] = new DoubleMatrix(Bc.viewRowSelection(blockGenes).getMatrix().toArray());
		}

		// Perform per block multiplication
		for (int i = 0; i < nrBlocks; i++) {
			// Re-cyle temorary results matrix
			DoubleMatrix Ctmp = new DoubleMatrix(Acache[i].rows, 1);
			int[] AcacheIndicesBlock = AcacheIndices[i];
			for (int col = 0; col < Bc.columns(); col++) {
				// Jblas multiplication
				Acache[i].mmuli(Bcache[i].getColumn(col), Ctmp);
				for (int j = 0; j < Ctmp.rows; j++) {
					int idx = AcacheIndicesBlock[j];
					C.setElementQuick(idx, col, C.getElementQuick(idx, col) + Ctmp.get(j, 0));
				}
				pb.step();
			}
		}
		pb.close();

		return C;
	}

	/**
	 * Calculate A %*% B for a matrix A whose structure has some form of block
	 * diagonality. I.e. there are rows and columns that are exactly zero.
	 *
	 * @param index gives the indices to multiply together.
	 * <p>
	 * Columns in A are first subset per block. The rows to keep in A are
	 * determined by removing rows that have only zero values in the subset of
	 * A's columns.
	 * <p>
	 * Indices provided in index should match the columns of A and the rows of
	 * B.
	 *
	 * @param Ac
	 * @param Bc
	 * @param index
	 * @return
	 */
	private static DoubleMatrixDataset<String, String> multBlockDiagonalSubsetPerColPcolt(final DoubleMatrixDataset<String, String> A, final DoubleMatrixDataset<String, String> B, final LinkedHashMap<String, ArrayList<String>> blockDiagonalIndices) {

		final DoubleMatrixDataset C = new DoubleMatrixDataset(A.getRowObjects(), B.getColObjects());

		final int nrBlocks = blockDiagonalIndices.size();

		ProgressBar pb = new ProgressBar("Multiplying per block", B.columns() * nrBlocks, ProgressBarStyle.ASCII);

		DoubleMatrix2D[] Acache = new DoubleMatrix2D[nrBlocks];
		DoubleMatrix2D[] Bcache = new DoubleMatrix2D[nrBlocks];

		int[][] AcacheIndices = new int[nrBlocks][];

		final ArrayList<String> blocks = new ArrayList(blockDiagonalIndices.keySet());

		// Compile the subset caches. Precomputing saves overhead. Views are returned so is memory efficient.
		for (int i = 0; i < nrBlocks; i++) {

			final ArrayList<String> blockGenes = blockDiagonalIndices.get(blocks.get(i));

			Acache[i] = A.viewColSelection(blockGenes).getMatrix();
			ArrayList<Integer> tmpIndices = new ArrayList<>();

			for (int j = 0; j < A.rows(); j++) {
				double rowSum = Acache[i].viewRow(j).zSum();
				if (rowSum != 0) {
					tmpIndices.add(j);
				}
			}

			AcacheIndices[i] = tmpIndices.stream().mapToInt(p -> p).toArray();
			Acache[i] = Acache[i].viewSelection(AcacheIndices[i], null);
			Bcache[i] = B.viewRowSelection(blockGenes).getMatrix();
		}

		// Perform per block multiplication
		for (int i = 0; i < nrBlocks; i++) {
			// Re-cyle temorary results matrix
			DoubleMatrix2D Ctmp = DoubleFactory2D.dense.make(Acache[i].rows(), 1);

			for (int col = 0; col < B.columns(); col++) {

				Acache[i].zMult(Bcache[i].viewSelection(null, new int[]{col}), Ctmp);
				for (int j = 0; j < Ctmp.rows(); j++) {
					int idx = AcacheIndices[i][j];
					C.setElementQuick(idx, col, C.getElementQuick(idx, col) + Ctmp.get(j, 0));
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
	 * Set mean to 0 and sd to 1 for each column in the data. Lifted from
	 * DoubleMatrixDataset.normalizeColumns() Refactored to center and scale to
	 * avoid confusion with DoubleMatrix2D.normalize() which normalizes the data
	 * so their sum = 1
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
	 * @param a Vector to calculate cumulative sum for
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
		Downstreamer.logDebugMem(msg, LOGGER);
	}

}
