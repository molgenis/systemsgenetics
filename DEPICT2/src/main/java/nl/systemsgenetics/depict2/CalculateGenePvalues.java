/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import cern.jet.math.tdouble.DoubleFunctions;
import com.opencsv.CSVWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.Executors;
import me.tongfei.progressbar.ProgressBar;
import me.tongfei.progressbar.ProgressBarStyle;
import static nl.systemsgenetics.depict2.Depict2.formatMsForLog;
import static nl.systemsgenetics.depict2.JamaHelperFunctions.eigenValueDecomposition;
import nl.systemsgenetics.depict2.originalLude.EstimateChi2SumDistUsingCorrelatedVariablesThread;
import org.apache.log4j.Logger;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetFastSubsetLoader;
import umcg.genetica.math.stats.ZScores;

/**
 *
 * @author patri
 */
public class CalculateGenePvalues {

	private static final Logger LOGGER = Logger.getLogger(CalculateGenePvalues.class);

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
	 * @param outputBasePath
	 * @return gene p-value matrix for each phenotype. rows: genes in same order
	 * as genes parameter, cols: phenotypes
	 * @throws java.io.IOException
	 */
	public static DoubleMatrixDataset<String, String> calculatorGenePvalues(
			final String variantPhenotypeZscoreMatrixPath,
			final GenotypeCorrelationSource genotypeCorrelationSource,
			final List<Gene> genes,
			final int windowExtend,
			final double maxR,
			final int nrPermutations,
			final String outputBasePath) throws IOException, Exception {

		final File geneVariantCountFile = new File(outputBasePath + "_geneVariantCount.txt");

		final List<String> phenotypes = Depict2.readMatrixAnnotations(new File(variantPhenotypeZscoreMatrixPath + ".cols.txt"));

		//Result matrix. Rows: genes, Cols: phenotypes
		final DoubleMatrixDataset<String, String> genePvalues = new DoubleMatrixDataset<>(createGeneHashRows(genes), createPhenoHashCols(phenotypes));
		final int numberPheno = phenotypes.size();
		final int numberGenes = genes.size();

		final DoubleMatrixDatasetFastSubsetLoader geneVariantPhenotypeMatrixRowLoader = new DoubleMatrixDatasetFastSubsetLoader(variantPhenotypeZscoreMatrixPath);

		final int[] genePValueDistributionPermuations = new int[21];//used to create histogram 
		final int[] genePValueDistributionChi2Dist = new int[21];//used to create histogram 

		int countRanPermutationsForGene = 0;
		int countBasedPvalueOnPermutations = 0;
		int countUseChi2DistForPvalue = 0;
		int countNoVariants = 0;

		long timeInPermutations = 0;
		long timeInCreatingGenotypeCorrelationMatrix = 0;
		long timeInDoingPca = 0;
		long timeInLoadingZscoreMatrix = 0;
		long timeInCalculatingPvalue = 0;

		long timeStart;
		long timeStop;

		final CSVWriter geneVariantCountWriter = new CSVWriter(new FileWriter(geneVariantCountFile), '\t', '\0', '\0', "\n");
		final String[] outputLine = new String[2];
		int c = 0;
		outputLine[c++] = "Gene";
		outputLine[c++] = "NumberVariants";
		geneVariantCountWriter.writeNext(outputLine);

		try (final ProgressBar pb = new ProgressBar("Gene p-value calculations", numberGenes, ProgressBarStyle.ASCII)) {

			for (int geneI = 0; geneI < numberGenes; ++geneI) {

				Gene gene = genes.get(geneI);

				timeStart = System.currentTimeMillis();

				final GenotypieCorrelationResult variantCorrelations = genotypeCorrelationSource.getCorrelationMatrixForRange(
						gene.getChr(), gene.getStart() - windowExtend, gene.getStop() + windowExtend, maxR);

				if (LOGGER.isDebugEnabled() & variantCorrelations.getIncludedVariants().length > 1) {

					DoubleMatrixDataset<String, String> corDataset = new DoubleMatrixDataset<>(variantCorrelations.getCorMatrix(), variantCorrelations.getIncludedVariants(), variantCorrelations.getIncludedVariants());

					corDataset.save(new File(outputBasePath + "_" + gene.getGene() + "_corMatrix.txt"));

				}

				timeStop = System.currentTimeMillis();
				timeInCreatingGenotypeCorrelationMatrix += (timeStop - timeStart);

				c = 0;
				outputLine[c++] = gene.getGene();
				outputLine[c++] = String.valueOf(variantCorrelations.getIncludedVariants().length);
				geneVariantCountWriter.writeNext(outputLine);

				final double[] geneChi2SumNull;
				if (variantCorrelations.getIncludedVariants().length > 1) {

					timeStart = System.currentTimeMillis();

					final Jama.EigenvalueDecomposition eig = eigenValueDecomposition(variantCorrelations.getCorMatrix());
					final double[] eigenValues = eig.getRealEigenvalues();

					if (LOGGER.isDebugEnabled()) {

						saveEigenValues(eigenValues, new File(outputBasePath + "_" + gene.getGene() + "_corMatrix.txt"));

					}

					timeStop = System.currentTimeMillis();
					timeInDoingPca += (timeStop - timeStart);

					timeStart = System.currentTimeMillis();

					geneChi2SumNull = runPermutationsUsingEigenValues(eigenValues, nrPermutations);

					timeStop = System.currentTimeMillis();
					timeInPermutations += (timeStop - timeStart);

					countRanPermutationsForGene++;
				} else {
					geneChi2SumNull = null;
				}

				timeStart = System.currentTimeMillis();

				//load current variants from variantPhenotypeMatrix
				final DoubleMatrixDataset<String, String> geneVariantPhenotypeMatrix = geneVariantPhenotypeMatrixRowLoader.loadSubsetOfRowsBinaryDoubleData(variantCorrelations.getIncludedVariants());

				timeStop = System.currentTimeMillis();
				timeInLoadingZscoreMatrix += (timeStop - timeStart);

				for (int phenoI = 0; phenoI < numberPheno; ++phenoI) {

					if (variantCorrelations.getIncludedVariants().length > 1) {

						timeStart = System.currentTimeMillis();

						final double geneChi2Sum = geneVariantPhenotypeMatrix.getCol(phenoI).aggregate(DoubleFunctions.plus, DoubleFunctions.square);

						double p = 0.5;
						for (int perm = 0; perm < nrPermutations; perm++) {
							if (geneChi2SumNull[perm] >= geneChi2Sum) {
								p++;
							}
						}
						p /= (double) nrPermutations + 1;

						if (p == 1) {
							p = 0.99999d;
						}
						if (p < 1e-300d) {
							p = 1e-300d;
						}

						genePValueDistributionPermuations[(int) (20d * p)]++;
						genePvalues.setElementQuick(geneI, phenoI, p);

						timeStop = System.currentTimeMillis();
						timeInCalculatingPvalue += (timeStop - timeStart);

						countBasedPvalueOnPermutations++;

					} else if (variantCorrelations.getIncludedVariants().length == 1) {

						//Always row 0
						double p = ZScores.zToP(-Math.abs(geneVariantPhenotypeMatrix.getElementQuick(0, phenoI)));
						if (p == 1) {
							p = 0.99999d;
						}
						if (p < 1e-300d) {
							p = 1e-300d;
						}
						genePValueDistributionChi2Dist[(int) (20d * p)]++;
						genePvalues.setElementQuick(geneI, phenoI, p);

						countUseChi2DistForPvalue++;

					} else {
						//no variants in or near gene
						//genePValueDistribution[(int) (20d * 0.99999d)]++;
						genePvalues.setElementQuick(geneI, phenoI, 0.99999d);
						countNoVariants++;
					}

				}

				pb.step();

			}

		}

		geneVariantCountWriter.close();

		LOGGER.debug("countRanPermutationsForGene: " + countRanPermutationsForGene);
		LOGGER.debug("countBasedPvalueOnPermutations: " + countBasedPvalueOnPermutations);
		LOGGER.debug("countUseChi2DistForPvalue: " + countUseChi2DistForPvalue);
		LOGGER.debug("countNoVariants: " + countNoVariants);

		LOGGER.debug("timeInPermutation: " + formatMsForLog(timeInPermutations));
		LOGGER.debug("timeInCreatingGenotypeCorrelationMatrix: " + formatMsForLog(timeInCreatingGenotypeCorrelationMatrix));
		LOGGER.debug("timeInDoingPca: " + formatMsForLog(timeInDoingPca));
		LOGGER.debug("timeInLoadingZscoreMatrix: " + formatMsForLog(timeInLoadingZscoreMatrix));
		LOGGER.debug("timeInCalculatingPvalue: " + formatMsForLog(timeInCalculatingPvalue));

		LOGGER.info("-----------------------");
		LOGGER.info("Gene p-value histrogram chi2 dist");
		for (double histCount : genePValueDistributionChi2Dist) {
			LOGGER.info(histCount);
		}
		LOGGER.info("Gene p-value histrogram permuations");
		for (double histCount : genePValueDistributionPermuations) {
			LOGGER.info(histCount);
		}
		LOGGER.info("-----------------------");

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

	private static void saveEigenValues(double[] eigenValues, File file) throws IOException {

		final CSVWriter eigenWriter = new CSVWriter(new FileWriter(file), '\t', '\0', '\0', "\n");
		final String[] outputLine = new String[2];
		int c = 0;
		outputLine[c++] = "Component";
		outputLine[c++] = "EigenValue";
		eigenWriter.writeNext(outputLine);

		for (int i = 0; i < eigenValues.length; ++i) {

			c = 0;
			outputLine[c++] = "PC" + i + 1;
			outputLine[c++] = String.valueOf(eigenValues[i]);
			eigenWriter.writeNext(outputLine);
			
		}

		eigenWriter.close();

	}

}
