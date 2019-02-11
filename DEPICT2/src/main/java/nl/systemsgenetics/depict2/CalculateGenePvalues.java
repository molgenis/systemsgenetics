/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import cern.jet.math.tdouble.DoubleFunctions;
import java.io.File;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.Executors;
import me.tongfei.progressbar.ProgressBar;
import me.tongfei.progressbar.ProgressBarStyle;
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
	 * @param variantPhenotypeZscoreMatrixPath Binary matrix. rows: variants,
	 * cols: phenotypes. Variants must have IDs identical to
	 * genotypeCovarianceSource
	 * @param genotypeCorrelationSource Source data used to calculate
	 * correlations between SNPs
	 * @param genes genes for which to calculate p-values.
	 * @param windowExtend number of bases to add left and right of gene window
	 * @param maxR max correlation between variants to use
	 * @param nrPermutations
	 * @return gene p-value matrix for each phenotype. rows: genes in same order
	 * as genes parameter, cols: phenotypes
	 * @throws java.io.IOException
	 */
	public static DoubleMatrixDataset<String, String> calculatorGenePvalues(
			final String variantPhenotypeZscoreMatrixPath,
			final GenotypeCovarianceSource genotypeCorrelationSource,
			final List<Gene> genes,
			final int windowExtend,
			final double maxR,
			final int nrPermutations) throws IOException {

		List<String> phenotypes = Depict2.readMatrixAnnotations(new File(variantPhenotypeZscoreMatrixPath + ".cols.txt"));

		//Result matrix. Rows: genes, Cols: phenotypes
		final DoubleMatrixDataset<String, String> genePvalues = new DoubleMatrixDataset<>(createGeneHashRows(genes), createPhenoHashCols(phenotypes));
		final int numberPheno = phenotypes.size();
		final int numberGenes = genes.size();

		final int[] genePValueDistribution = new int[21];//used to create histogram 

		try (ProgressBar pb = new ProgressBar("Gene p-value calculations", 10, ProgressBarStyle.ASCII)) {

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

				//load current variants from variantPhenotypeMatrix
				final DoubleMatrixDataset<String, String> geneVariantPhenotypeMatrix = DoubleMatrixDataset.loadSubsetOfRowsBinaryDoubleData(variantPhenotypeZscoreMatrixPath, variantCorrelations.getIncludedVariants());

				for (int phenoI = 0; phenoI < numberPheno; ++phenoI) {

					if (variantCorrelations.getIncludedVariants().length > 1) {

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

						//Always row 0
						double p = ZScores.zToP(-Math.abs(geneVariantPhenotypeMatrix.getElementQuick(0, phenoI)));
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

				pb.step();

//			if (geneI % 100 == 0 & geneI > 0) {
//				System.out.print("Proccessed " + geneI + " genes");
//			}
			}

		}

		System.out.println("-----------------------");
		System.out.println("Gene p-value histrogram");
		for (double histCount : genePValueDistribution) {
			System.out.println(histCount);
		}
		System.out.println("-----------------------");

		return genePvalues;

	}

	public static LinkedHashMap<String, Integer> createGeneHashRows(final List<Gene> genes) {
		LinkedHashMap<String, Integer> geneHashRows = new LinkedHashMap<>(genes.size());
		for (int geneI = 0; geneI < genes.size(); ++geneI) {
			geneHashRows.put(genes.get(geneI).getGene(), geneI);
		}
		return geneHashRows;
	}

	public static LinkedHashMap<String, Integer> createPhenoHashCols(final List<String> phenotypes) {
		LinkedHashMap<String, Integer> phenoHashCols = new LinkedHashMap<>(phenotypes.size());
		for (int phenoI = 0; phenoI < phenotypes.size(); ++phenoI) {
			phenoHashCols.put(phenotypes.get(phenoI), phenoI);
		}
		return phenoHashCols;
	}

	public static double[] runPermutationsUsingEigenValues(final double[] eigenValues, final int nrPermutations) {

		final int nrTasks = 1000;
		final int nrPermutationsPerTask = nrPermutations / nrTasks;
		final double[] geneChi2SumNull = new double[nrPermutations];

		int nrThreads = Depict2Options.getNumberOfThreadsToUse();
		java.util.concurrent.ExecutorService threadPool = Executors.newFixedThreadPool(nrThreads);
		CompletionService<double[]> pool = new ExecutorCompletionService<>(threadPool);
		for (int taskNr = 0; taskNr < nrTasks; taskNr++) {
			EstimateChi2SumDistUsingCorrelatedVariablesThread task = new EstimateChi2SumDistUsingCorrelatedVariablesThread(eigenValues, nrPermutationsPerTask, taskNr);
			pool.submit(task);
		}
		try {
			for (int task = 0; task < nrTasks; task++) {
				try {
					double[] chi2SumDist = pool.take().get();
					System.arraycopy(chi2SumDist, 0, geneChi2SumNull, task * nrPermutationsPerTask, chi2SumDist.length);
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
