/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import cern.jet.math.tdouble.DoubleFunctions;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.Executors;
import static nl.systemsgenetics.depict2.JamaHelperFunctions.eigenValueDecomposition;
import nl.systemsgenetics.depict2.originalLude.DoubleArrayIntegerObject;
import nl.systemsgenetics.depict2.originalLude.EstimateChi2SumDistUsingCorrelatedVariablesThread;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.ZScores;

/**
 *
 * @author patri
 */
public class CalculateGenePvalues {

	/**
	 *
	 * @param variantPhenotypeZscoreMatrix rows: variants, cols: phenotypes.
	 * Variants must have IDs identical to genotypeCovarianceSource
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
			final DoubleMatrixDataset<String, String> variantPhenotypeZscoreMatrix,
			final GenotypeCovarianceSource genotypeCorrelationSource,
			final ArrayList<Gene> genes,
			final int windowExtend,
			final double maxR,
			final int nrPermutations) {

		//Result matrix. Rows: genes, Cols: phenotypes
		final DoubleMatrixDataset<String, String> genePvalues = new DoubleMatrixDataset<>(createGeneHashRows(genes), variantPhenotypeZscoreMatrix.getHashColsCopy());
		final int numberPheno = variantPhenotypeZscoreMatrix.columns();
		final int numberGenes = variantPhenotypeZscoreMatrix.rows();

		final int[] genePValueDistribution = new int[21];//used to create histogram 

		for (int geneI = 0; geneI < numberGenes; ++geneI) {

			Gene gene = genes.get(geneI);

			final GenotypieCorrelationResult variantCorrelations = genotypeCorrelationSource.getCorrelationMatrixForRange(
					gene.getChr(), gene.getStart() - windowExtend, gene.getStop() + windowExtend, maxR);

			final double[] geneChi2SumNull;

			if (variantCorrelations.getIncludedVariants().length > 1) {
				final Jama.EigenvalueDecomposition eig = eigenValueDecomposition(variantCorrelations.getCorMatrix());
				final double[] eigenValues = eig.getRealEigenvalues();
				geneChi2SumNull = runPermutationsUsingEigenValues(eigenValues, nrPermutations);
			} else {
				geneChi2SumNull = null;
			}

			for (int phenoI = 0; phenoI < numberPheno; ++phenoI) {

				if (variantCorrelations.getIncludedVariants().length > 1) {

					//extract current variants from variantPhenotypeMatrix
					final DoubleMatrixDataset<String, String> geneVariantPhenotypeMatrix = variantPhenotypeZscoreMatrix.viewRowSelection(variantCorrelations.getIncludedVariants());

					final double geneChi2Sum = geneVariantPhenotypeMatrix.getCol(phenoI).aggregate(DoubleFunctions.plus, DoubleFunctions.square);

					double p = 0.5;
					for (int perm = 0; perm < nrPermutations; perm++) {
						if (geneChi2SumNull[perm] >= geneChi2Sum) {
							p += 1;
						}
					}
					p /= (double) nrPermutations + 1;
					
					if (p == 1) {
						p = 0.99999d;
					}
					if (p < 1e-300d) {
						p = 1e-300d;
					}
					
					genePValueDistribution[(int) (20d * p)]++;
					genePvalues.setElementQuick(geneI, phenoI, p);

				} else if (variantCorrelations.getIncludedVariants().length == 1) {

					final int variantI = variantPhenotypeZscoreMatrix.getRowIndex(variantCorrelations.getIncludedVariants()[0]);

					double p = ZScores.zToP(-Math.abs(variantPhenotypeZscoreMatrix.getElementQuick(variantI, phenoI)));
					if (p == 1) {
						p = 0.99999d;
					}
					if (p < 1e-300d) {
						p = 1e-300d;
					}
					genePValueDistribution[(int) (20d * p)]++;
					genePvalues.setElementQuick(geneI, phenoI, p);

				} else {
					//no variants in or near gene
					genePValueDistribution[(int) (20d * 0.99999d)]++;
					genePvalues.setElementQuick(geneI, phenoI, 0.99999d);
				}

			}

		}

		System.out.println("-----------------------");
		System.out.println("Gene p-value histrogram");
		for(double histCount : genePValueDistribution){
			System.out.println(histCount);
		}
		System.out.println("-----------------------");
		
		return genePvalues;

	}

	public static LinkedHashMap<String, Integer> createGeneHashRows(final ArrayList<Gene> genes) {
		LinkedHashMap<String, Integer> geneHashRows = new LinkedHashMap<>(genes.size());
		for (int geneI = 0; geneI < genes.size(); ++geneI) {
			geneHashRows.put(genes.get(geneI).getGene(), geneI);
		}
		return geneHashRows;
	}

	public static double[] runPermutationsUsingEigenValues(final double[] eigenValues, final int nrPermutations) {

		final int nrTasks = 1000;
		final int nrPermutationsPerTask = nrPermutations / nrTasks;
		final double[] geneChi2SumNull = new double[nrPermutations];

		int nrThreads = Depict2Options.getNumberOfThreadsToUse();
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
					System.arraycopy(chi2SumDist, 0, geneChi2SumNull, taskNr * nrPermutationsPerTask, chi2SumDist.length);
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
