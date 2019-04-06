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
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.Level;
import java.util.stream.IntStream;
import me.tongfei.progressbar.ProgressBar;
import me.tongfei.progressbar.ProgressBarStyle;
import static nl.systemsgenetics.depict2.Depict2.formatMsForLog;
import static nl.systemsgenetics.depict2.JamaHelperFunctions.eigenValueDecomposition;
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
	protected static long timeInCreatingGenotypeCorrelationMatrix = 0;
	protected static long timeInLoadingGenotypeDosages = 0;

	static long timeInPermutations = 0;
	static long timeInPruningGenotypeCorrelationMatrix = 0;
	static long timeInDoingPca = 0;
	static long timeInLoadingZscoreMatrix = 0;
	static long timeInCalculatingPvalue = 0;
	static long timeInCalculatingRealSumChi2 = 0;

	static int countRanPermutationsForGene = 0;
	static int countBasedPvalueOnPermutations = 0;
	static int countUseChi2DistForPvalue = 0;
	static int countNoVariants = 0;

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
	 * @param randomChi2
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
			final String outputBasePath,
			final double[] randomChi2) throws IOException, Exception {

		final File geneVariantCountFile = new File(outputBasePath + "_geneVariantCount.txt");

		final List<String> phenotypes = Depict2.readMatrixAnnotations(new File(variantPhenotypeZscoreMatrixPath + ".cols.txt"));

		//Result matrix. Rows: genes, Cols: phenotypes
		final DoubleMatrixDataset<String, String> genePvalues = new DoubleMatrixDataset<>(createGeneHashRows(genes), createPhenoHashCols(phenotypes));
		final int numberPheno = phenotypes.size();
		final int numberGenes = genes.size();

		final DoubleMatrixDatasetFastSubsetLoader geneVariantPhenotypeMatrixRowLoader = new DoubleMatrixDatasetFastSubsetLoader(variantPhenotypeZscoreMatrixPath);

		final int[] genePValueDistributionPermuations = new int[21];//used to create histogram 
		final int[] genePValueDistributionChi2Dist = new int[21];//used to create histogram 

		final CSVWriter geneVariantCountWriter = new CSVWriter(new FileWriter(geneVariantCountFile), '\t', '\0', '\0', "\n");
		final String[] outputLine = new String[2];
		int c = 0;
		outputLine[c++] = "Gene";
		outputLine[c++] = "NumberVariants";
		geneVariantCountWriter.writeNext(outputLine);

		try (final ProgressBar pb = new ProgressBar("Gene p-value calculations", numberGenes, ProgressBarStyle.ASCII)) {

			//All genes are indipe
			IntStream.range(0, numberGenes).parallel().forEach((int geneI) -> {

				try {
					runGene(genes, geneI, genotypeCorrelationSource, windowExtend, maxR, outputBasePath, geneVariantCountWriter, nrPermutations, randomChi2, geneVariantPhenotypeMatrixRowLoader, numberPheno, genePValueDistributionPermuations, genePvalues, genePValueDistributionChi2Dist);
					pb.step();
				} catch (Exception ex) {
					throw new RuntimeException(ex);
				}

			});

		}

		geneVariantCountWriter.close();

		LOGGER.debug("countRanPermutationsForGene: " + countRanPermutationsForGene);
		LOGGER.debug("countBasedPvalueOnPermutations: " + countBasedPvalueOnPermutations);
		LOGGER.debug("countUseChi2DistForPvalue: " + countUseChi2DistForPvalue);
		LOGGER.debug("countNoVariants: " + countNoVariants);

		LOGGER.info("timeInLoadingGenotypeDosages: " + formatMsForLog(timeInLoadingGenotypeDosages));
		LOGGER.info("timeInCreatingGenotypeCorrelationMatrix: " + formatMsForLog(timeInCreatingGenotypeCorrelationMatrix));
		LOGGER.info("timeInPruningGenotypeCorrelationMatrix: " + formatMsForLog(timeInPruningGenotypeCorrelationMatrix));
		LOGGER.info("timeInDoingPca: " + formatMsForLog(timeInDoingPca));
		LOGGER.info("timeInPermutation: " + formatMsForLog(timeInPermutations));
		LOGGER.info("timeInLoadingZscoreMatrix: " + formatMsForLog(timeInLoadingZscoreMatrix));
		LOGGER.info("timeInCalculatingRealSumChi2: " + formatMsForLog(timeInCalculatingRealSumChi2));
		LOGGER.info("timeInCalculatingPvalue: " + formatMsForLog(timeInCalculatingPvalue));

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

	private static void runGene(final List<Gene> genes, int geneI, final GenotypeCorrelationSource genotypeCorrelationSource, final int windowExtend, final double maxR, final String outputBasePath, final CSVWriter geneVariantCountWriter, final int nrPermutations, final double[] randomChi2, final DoubleMatrixDatasetFastSubsetLoader geneVariantPhenotypeMatrixRowLoader, final int numberPheno, final int[] genePValueDistributionPermuations, final DoubleMatrixDataset<String, String> genePvalues, final int[] genePValueDistributionChi2Dist) throws Exception, IOException {

		long timeStart;
		long timeStop;

		Gene gene = genes.get(geneI);

		final DoubleMatrixDataset<String, String> variantCorrelationsPruned;

		{
			//timeStart = System.currentTimeMillis();
			final DoubleMatrixDataset<String, String> variantCorrelations = genotypeCorrelationSource.getCorrelationMatrixForRange(
					gene.getChr(), gene.getStart() - windowExtend, gene.getStop() + windowExtend);
//					timeStop = System.currentTimeMillis();
//					timeInCreatingGenotypeCorrelationMatrix += (timeStop - timeStart);

			timeStart = System.currentTimeMillis();
			variantCorrelationsPruned = pruneCorrelatinMatrix(variantCorrelations, maxR);
			timeStop = System.currentTimeMillis();
			timeInPruningGenotypeCorrelationMatrix += (timeStop - timeStart);
		}

		if (LOGGER.isDebugEnabled() & variantCorrelationsPruned.rows() > 1) {

			variantCorrelationsPruned.save(new File(outputBasePath + "_" + gene.getGene() + "_corMatrix.txt"));

		}

		timeStop = System.currentTimeMillis();
		timeInCreatingGenotypeCorrelationMatrix += (timeStop - timeStart);

		int c = 0;
		final String[] outputLine = new String[2];
		outputLine[c++] = gene.getGene();
		outputLine[c++] = String.valueOf(variantCorrelationsPruned.rows());
		//This is writeNext is threadsafe
		geneVariantCountWriter.writeNext(outputLine);

		final double[] geneChi2SumNull;
		if (variantCorrelationsPruned.rows() > 1) {

			timeStart = System.currentTimeMillis();

			final Jama.EigenvalueDecomposition eig = eigenValueDecomposition(variantCorrelationsPruned.getMatrixAs2dDoubleArray());
			final double[] eigenValues = eig.getRealEigenvalues();
			for (int e = 0; e < eigenValues.length; e++) {
				if (eigenValues[e] < 0) {
					eigenValues[e] = 0;
				}
			}

			if (LOGGER.isDebugEnabled()) {

				saveEigenValues(eigenValues, new File(outputBasePath + "_" + gene.getGene() + "_eigenValues.txt"));

			}

			timeStop = System.currentTimeMillis();
			timeInDoingPca += (timeStop - timeStart);

			timeStart = System.currentTimeMillis();

			geneChi2SumNull = runPermutationsUsingEigenValues(eigenValues, nrPermutations, randomChi2);

			timeStop = System.currentTimeMillis();
			timeInPermutations += (timeStop - timeStart);

			countRanPermutationsForGene++;
		} else {
			geneChi2SumNull = null;
		}

		timeStart = System.currentTimeMillis();

//load current variants from variantPhenotypeMatrix
		final DoubleMatrixDataset<String, String> geneVariantPhenotypeMatrix = geneVariantPhenotypeMatrixRowLoader.loadSubsetOfRowsBinaryDoubleData(variantCorrelationsPruned.getRowObjects());

		timeStop = System.currentTimeMillis();
		timeInLoadingZscoreMatrix += (timeStop - timeStart);

		for (int phenoI = 0; phenoI < numberPheno; ++phenoI) {

			if (variantCorrelationsPruned.rows() > 1) {

				timeStart = System.currentTimeMillis();

				final double geneChi2Sum = geneVariantPhenotypeMatrix.getCol(phenoI).aggregate(DoubleFunctions.plus, DoubleFunctions.square);

				timeStop = System.currentTimeMillis();
				timeInCalculatingRealSumChi2 += (timeStop - timeStart);

				timeStart = System.currentTimeMillis();

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

			} else if (variantCorrelationsPruned.rows() == 1) {

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
				genePvalues.setElementQuick(geneI, phenoI, 0.5d);
				countNoVariants++;
			}

		}

		
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

	public static double[] runPermutationsUsingEigenValues(final double[] eigenValues, final int nrPermutations, double[] randomChi2) {

		final double[] geneChi2SumNull = new double[nrPermutations];

		final int randomChi2Size = randomChi2.length;
		ThreadLocalRandom rnd = ThreadLocalRandom.current();

		int i = 0;
		for (int perm = 0; perm < nrPermutations; perm++) {
			double weightedChi2Perm = 0;
			for (int g = 0; g < eigenValues.length; g++) {

				if (i < randomChi2Size) {
					weightedChi2Perm += eigenValues[g] * randomChi2[i++];
				} else {
					double randomZ = rnd.nextGaussian();
					weightedChi2Perm += eigenValues[g] * randomZ * randomZ;
				}

			}
			geneChi2SumNull[perm] = weightedChi2Perm;
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
			outputLine[c++] = "PC" + (i + 1);
			outputLine[c++] = String.valueOf(eigenValues[i]);
			eigenWriter.writeNext(outputLine);

		}

		eigenWriter.close();

	}

	private static DoubleMatrixDataset<String, String> pruneCorrelatinMatrix(DoubleMatrixDataset<String, String> correlationMatrixForRange, double maxR) {

		ArrayList<String> variantNames = correlationMatrixForRange.getRowObjects();
		LinkedHashSet<String> includedVariants = new LinkedHashSet<>(correlationMatrixForRange.rows());

		rows:
		for (int r = 0; r < correlationMatrixForRange.rows(); ++r) {
			cols:
			for (int c = 0; c < r; ++c) {
				if (Math.abs(correlationMatrixForRange.getElementQuick(r, c)) >= maxR && includedVariants.contains(variantNames.get(c))) {
					continue rows;
				}
			}
			includedVariants.add(variantNames.get(r));
		}

		LOGGER.debug(" * Variants after pruning high r: " + includedVariants.size());

		return correlationMatrixForRange.viewSelection(includedVariants, includedVariants);

	}

}
