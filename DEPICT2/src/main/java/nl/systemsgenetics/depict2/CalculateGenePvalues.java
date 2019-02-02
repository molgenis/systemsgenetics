/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import java.util.ArrayList;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.Executors;
import static nl.systemsgenetics.depict2.JamaHelperFunctions.eigenValueDecomposition;
import nl.systemsgenetics.depict2.originalLude.DoubleArrayIntegerObject;
import nl.systemsgenetics.depict2.originalLude.EstimateChi2SumDistUsingCorrelatedVariablesThread;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class CalculateGenePvalues {

	/**
	 *
	 * @param variantPhenotypeMatrix rows: variants, cols: phenotypes. Variants
	 * must have IDs identical to genotypeCovarianceSource
	 * @param genotypeCorrelationSource Source data used to calculate
	 * correlations between SNPs
	 * @param genes genes for which to calculate p-values.
	 * @param windowExtend number of bases to add left and right of gene window
	 * @param maxR max correlation between variants to use
	 * @param nrPermutations 
	 * @return gene p-value matrix for each phenotype. rows: genes in same order
	 * as genes parameter, cols: phenotypes
	 */
	public static DoubleMatrixDataset<String, String> calculatorGenePvalues(
			final DoubleMatrixDataset<String, String> variantPhenotypeMatrix,
			final GenotypeCovarianceSource genotypeCorrelationSource,
			final ArrayList<Gene> genes,
			final int windowExtend,
			final double maxR,
			final int nrPermutations) {

		for (Gene gene : genes) {

			final GenotypieCorrelationResult variantCorrelations = genotypeCorrelationSource.getCorrelationMatrixForRange(
					gene.getChr(), gene.getStart() - windowExtend, gene.getStop() + windowExtend, maxR);

			Jama.EigenvalueDecomposition eig = eigenValueDecomposition(variantCorrelations.getCorMatrix());
			double[] eigenValues = eig.getRealEigenvalues();
			
			final double[] geneChi2SumNull = runPermutationsUsingEigenValues(eigenValues, nrPermutations);

			//extract current variants from variantPhenotypeMatrix
			final DoubleMatrixDataset<String, String> geneVariantPhenotypeMatrix = variantPhenotypeMatrix.viewRowSelection(variantCorrelations.getIncludedVariants());

		}

		return null;

	}

	public static double[] runPermutationsUsingEigenValues(final double[] eigenValues, final int nrPermutations) {

		final int nrTasks = 1000;
		final int nrPermutationsPerTask = nrPermutations / nrTasks;
		final double[] geneChi2SumNull = new double[nrPermutations];

		int nrThreads = Runtime.getRuntime().availableProcessors();
		java.util.concurrent.ExecutorService threadPool = Executors.newFixedThreadPool(nrThreads);
		CompletionService<DoubleArrayIntegerObject> pool = new ExecutorCompletionService<>(threadPool);
		for (int taskNr = 0; taskNr < nrTasks; taskNr++) {
			EstimateChi2SumDistUsingCorrelatedVariablesThread task = new EstimateChi2SumDistUsingCorrelatedVariablesThread(eigenValues, nrPermutationsPerTask, taskNr);
			pool.submit(task);
		}
		try {
			for (int task = 0; task < nrTasks; task++) {
				try {
					DoubleArrayIntegerObject result = pool.take().get();
					int taskNr = result.intValue;
					double[] chi2SumDist = result.doubleArray;
					System.arraycopy(chi2SumDist, 0, geneChi2SumNull, taskNr * 10000, chi2SumDist.length);
				} catch (ExecutionException ex) {
					throw new RuntimeException(ex);
				}
			}
			threadPool.shutdown();
		} catch (InterruptedException | RuntimeException e) {
			throw new RuntimeException(e);
		}

		return geneChi2SumNull;

	}

}
