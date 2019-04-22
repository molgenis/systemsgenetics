/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.jet.math.tdouble.DoubleFunctions;
import com.opencsv.CSVWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.concurrent.ThreadLocalRandom;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.logging.Level;
import java.util.stream.IntStream;
import me.tongfei.progressbar.ProgressBar;
import me.tongfei.progressbar.ProgressBarStyle;
import static nl.systemsgenetics.depict2.Depict2.formatMsForLog;
import static nl.systemsgenetics.depict2.JamaHelperFunctions.eigenValueDecomposition;
import org.apache.log4j.Logger;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetFastSubsetLoader;
import umcg.genetica.math.stats.PearsonRToZscoreBinned;
import umcg.genetica.math.stats.ZScores;

/**
 *
 * @author patri
 */
public class GenePvalueCalculator {

	private static final Logger LOGGER = Logger.getLogger(GenePvalueCalculator.class);
	private static final DoubleMatrixDataset<String, String> EMPTY_DATASET = new DoubleMatrixDataset<>(0, 0);
	private static final int PERMUTATION_STEP = 10000;
	private static final int NUMBER_RANDOM_PHENO = PERMUTATION_STEP;
	private static final int NUMBER_PERMUTATION_NULL_GWAS = 10000;
	private static final double NUMBER_PERMUTATION_NULL_GWAS_PLUS_1 = NUMBER_PERMUTATION_NULL_GWAS + 1;

	protected static long timeInCreatingGenotypeCorrelationMatrix = 0;
	protected static long timeInLoadingGenotypeDosages = 0;

	private static long totalTimeInThread = 0;
	private static long timeInPermutations = 0;
	private static long timeInPruningGenotypeCorrelationMatrix = 0;
	private static long timeInDoingPca = 0;
	private static long timeInLoadingZscoreMatrix = 0;
	private static long timeInCalculatingPvalue = 0;
	private static long timeInCalculatingRealSumChi2 = 0;
	private static long timeInComparingRealChi2ToPermutationChi2 = 0;
	private static long timeInNormalizingGenotypeDosages = 0;

	private static int countRanPermutationsForGene = 0;
	private static int countBasedPvalueOnPermutations = 0;
	private static int countUseChi2DistForPvalue = 0;
	private static int countNoVariants = 0;

	private final RandomAccessGenotypeData referenceGenotypes;
	private final List<Gene> genes;
	private final int windowExtend;
	private final double maxR;
	private final int maxNrPermutations;
	private final double nrPermutationsPlus1Double;
	private final String outputBasePath;
	private final double[] randomChi2;
	private final LinkedHashMap<String, Integer> sampleHash;//In order of the samples in the genotype data
	private final DoubleMatrixDataset<String, String> randomNormalizedPhenotypes;//Rows samples order of sample hash, cols phenotypes
	private final PearsonRToZscoreBinned r2zScore;//Based on number of samples with genotypes
	private final DoubleMatrixDataset<String, String> genePvalues;
	private final DoubleMatrixDataset<String, String> genePvaluesNullGwas;
	private final DoubleMatrixDataset<String, String> geneVariantCount;
	private final DoubleMatrixDataset<String, String> geneMaxPermutationCount;
	private final int numberRealPheno;
	private final DoubleMatrixDatasetFastSubsetLoader geneVariantPhenotypeMatrixRowLoader;
	private final int[] genePValueDistributionPermuations;
	private final int[] genePValueDistributionChi2Dist;
	private final ProgressBar pb;

	/**
	 *
	 * @param variantPhenotypeZscoreMatrixPath Binary matrix. rows: variants,
	 * cols: phenotypes. Variants must have IDs identical to
	 * genotypeCovarianceSource
	 * @param referenceGenotypes Source data used to calculate correlations
	 * between SNPs
	 * @param genes genes for which to calculate p-values.
	 * @param windowExtend number of bases to add left and right of gene window
	 * @param maxR max correlation between variants to use
	 * @param nrPermutations
	 * @param outputBasePath
	 * @param randomChi2
	 * @throws java.io.IOException
	 */
	@SuppressWarnings("CallToThreadStartDuringObjectConstruction")
	public GenePvalueCalculator(String variantPhenotypeZscoreMatrixPath, RandomAccessGenotypeData referenceGenotypes, List<Gene> genes, int windowExtend, double maxR, int nrPermutations, String outputBasePath, double[] randomChi2) throws IOException, Exception {

		this.referenceGenotypes = referenceGenotypes;
		this.genes = genes;
		this.windowExtend = windowExtend;
		this.maxR = maxR;
		this.maxNrPermutations = nrPermutations;
		this.nrPermutationsPlus1Double = nrPermutations + 1;
		this.outputBasePath = outputBasePath;
		this.randomChi2 = randomChi2;

		if (nrPermutations % PERMUTATION_STEP != 0) {
			throw new Exception("Number of permutaitons must be dividable by: " + PERMUTATION_STEP);
		}

		String[] samples = referenceGenotypes.getSampleNames();

		sampleHash = new LinkedHashMap<>(samples.length);

		int s = 0;
		for (String sample : samples) {
			sampleHash.put(sample, s++);
		}

		randomNormalizedPhenotypes = generateRandomNormalizedPheno(sampleHash);

		r2zScore = new PearsonRToZscoreBinned(10000000, sampleHash.size());//10000000

		final List<String> phenotypes = Depict2.readMatrixAnnotations(new File(variantPhenotypeZscoreMatrixPath + ".cols.txt"));

		//Result matrix. Rows: genes, Cols: phenotypes
		genePvalues = new DoubleMatrixDataset<>(createGeneHashRows(genes), createHashColsFromList(phenotypes));

		//Result matrix null GWAS. Rows: genes, Cols: phenotypes
		genePvaluesNullGwas = new DoubleMatrixDataset<>(createGeneHashRows(genes), randomNormalizedPhenotypes.getHashCols());

		geneVariantCount = new DoubleMatrixDataset<>(createGeneHashRows(genes), createHashColsFromList(Arrays.asList(new String[]{"count"})));
		geneMaxPermutationCount = new DoubleMatrixDataset<>(createGeneHashRows(genes), createHashColsFromList(Arrays.asList(new String[]{"maxPermutations"})));

		numberRealPheno = phenotypes.size();
		final int numberGenes = genes.size();

		geneVariantPhenotypeMatrixRowLoader = new DoubleMatrixDatasetFastSubsetLoader(variantPhenotypeZscoreMatrixPath);

		genePValueDistributionPermuations = new int[21];//used to create histogram 
		genePValueDistributionChi2Dist = new int[21];//used to create histogram 

		pb = new ProgressBar("Gene p-value calculations", numberGenes, ProgressBarStyle.ASCII);

		//All genes are indipendant
		// IntStream.range(0, numberGenes).parallel().forEach((int geneI) -> {
//			for (int geneI = 0; geneI < numberGenes; ++geneI) {
		long startThread = System.currentTimeMillis();

		//Results are writen in genePvalues
		//runGene(geneI);
		final AtomicInteger count = new AtomicInteger(0);
		List<Thread> threads = new ArrayList<>(Depict2Options.getNumberOfThreadsToUse());
		Depict2.ThreadErrorHandler threadErrorHandler = new Depict2.ThreadErrorHandler("gene p-value calculator");

		for (int i = 0; i < Depict2Options.getNumberOfThreadsToUse(); ++i) {

			Thread worker = new Thread(new calculatorThread(count));
			worker.start();
			worker.setUncaughtExceptionHandler(threadErrorHandler);
			threads.add(worker);

		}

		boolean running;
		do {
			running = false;
			for (Thread thread : threads) {
				if (thread.isAlive()) {
					running = true;
				}
			}

			try {
				Thread.sleep(500);
			} catch (InterruptedException ex) {
			}

		} while (running);

		long endThread = System.currentTimeMillis();

		totalTimeInThread += (endThread - startThread);

		pb.close();

		LOGGER.debug("countRanPermutationsForGene: " + countRanPermutationsForGene);
		LOGGER.debug("countBasedPvalueOnPermutations: " + countBasedPvalueOnPermutations);
		LOGGER.debug("countUseChi2DistForPvalue: " + countUseChi2DistForPvalue);
		LOGGER.debug("countNoVariants: " + countNoVariants);

		LOGGER.info("timeInLoadingGenotypeDosages: " + formatMsForLog(timeInLoadingGenotypeDosages));
		LOGGER.info("timeInNormalizingGenotypeDosages: " + formatMsForLog(timeInNormalizingGenotypeDosages));
		LOGGER.info("timeInCreatingGenotypeCorrelationMatrix: " + formatMsForLog(timeInCreatingGenotypeCorrelationMatrix));
		LOGGER.info("timeInPruningGenotypeCorrelationMatrix: " + formatMsForLog(timeInPruningGenotypeCorrelationMatrix));
		LOGGER.info("timeInDoingPca: " + formatMsForLog(timeInDoingPca));
		LOGGER.info("timeInPermutation: " + formatMsForLog(timeInPermutations));
		LOGGER.info("timeInLoadingZscoreMatrix: " + formatMsForLog(timeInLoadingZscoreMatrix));
		LOGGER.info("timeInCalculatingRealSumChi2: " + formatMsForLog(timeInCalculatingRealSumChi2));
		LOGGER.info("timeInComparingRealChi2ToPermutationChi2: " + formatMsForLog(timeInComparingRealChi2ToPermutationChi2));
		LOGGER.info("timeInCalculatingPvalue: " + formatMsForLog(timeInCalculatingPvalue));
		LOGGER.info("totalTimeInThread: " + formatMsForLog(totalTimeInThread) + " (" + totalTimeInThread + ") ");

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

	}

	/**
	 *
	 * @return gene p-value matrix for each phenotype. rows: genes in same order
	 * as genes list, cols: phenotypes
	 * @throws java.io.IOException
	 */
	public DoubleMatrixDataset<String, String> getGenePvalues() {
		return genePvalues;
	}

	public DoubleMatrixDataset<String, String> getGenePvaluesNullGwas(){
		return genePvaluesNullGwas;
	}

	public DoubleMatrixDataset<String, String> getGeneVariantCount() {
		return geneVariantCount;
	}

	public DoubleMatrixDataset<String, String> getGeneMaxPermutationCount() {
		return geneMaxPermutationCount;
	}

	private void runGene(int geneI, final double[] geneChi2SumNull) throws IOException, Exception {
		long timeStart;
		long timeStop;

		Gene gene = genes.get(geneI);

		//Rows: samples, cols: variants
		final DoubleMatrixDataset<String, String> variantScaledDosagesPruned;
		final DoubleMatrixDataset<String, String> variantCorrelationsPruned;

		{

			DoubleMatrixDataset<String, String> variantScaledDosages = loadVariantScaledDosageMatrix(gene.getChr(), gene.getStart() - windowExtend, gene.getStop() + windowExtend);

			timeStart = System.currentTimeMillis();
			final DoubleMatrixDataset<String, String> variantCorrelations = variantScaledDosages.calculateCorrelationMatrixOnNormalizedColumns();
			timeStop = System.currentTimeMillis();
			timeInCreatingGenotypeCorrelationMatrix += (timeStop - timeStart);

			timeStart = System.currentTimeMillis();
			variantCorrelationsPruned = pruneCorrelatinMatrix(variantCorrelations, maxR);
			timeStop = System.currentTimeMillis();
			timeInPruningGenotypeCorrelationMatrix += (timeStop - timeStart);

			variantScaledDosagesPruned = variantScaledDosages.viewColSelection(variantCorrelationsPruned.getHashCols().keySet());

		}

		final DoubleMatrixDataset<String, String> nullGwasZscores;
		if (variantCorrelationsPruned.rows() > 0) {
			//nullGwasZscores will first contain peason r valus but this will be converted to Z-score using a lookup table;
			nullGwasZscores = DoubleMatrixDataset.correlateColumnsOf2ColumnNormalizedDatasets(randomNormalizedPhenotypes, variantScaledDosagesPruned);
			r2zScore.inplaceRToZ(nullGwasZscores);
		} else {
			nullGwasZscores = EMPTY_DATASET;
		}

		if (LOGGER.isDebugEnabled() & variantCorrelationsPruned.rows() > 1) {

			variantCorrelationsPruned.save(new File(outputBasePath + "_" + gene.getGene() + "_corMatrix.txt"));

		}

		timeStop = System.currentTimeMillis();
		timeInCreatingGenotypeCorrelationMatrix += (timeStop - timeStart);

		geneVariantCount.setElementQuick(geneI, 0, variantCorrelationsPruned.rows());

		int currentNumberPermutationsCalculated = 0;
		final ThreadLocalRandom rnd = ThreadLocalRandom.current();
		final double[] eigenValues;
		if (variantCorrelationsPruned.rows() > 1) {

			timeStart = System.currentTimeMillis();

			final Jama.EigenvalueDecomposition eig = eigenValueDecomposition(variantCorrelationsPruned.getMatrixAs2dDoubleArray());
			eigenValues = eig.getRealEigenvalues();
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

			timeStop = System.currentTimeMillis();
			timeInPermutations += (timeStop - timeStart);

			countRanPermutationsForGene++;
		} else {
			eigenValues = null;
		}

		timeStart = System.currentTimeMillis();

		//load current variants from variantPhenotypeMatrix
		final DoubleMatrixDataset<String, String> geneVariantPhenotypeMatrix;
		synchronized (this) {
			geneVariantPhenotypeMatrix = geneVariantPhenotypeMatrixRowLoader.loadSubsetOfRowsBinaryDoubleData(variantCorrelationsPruned.getRowObjects());
		}

		timeStop = System.currentTimeMillis();
		timeInLoadingZscoreMatrix += (timeStop - timeStart);

		for (int phenoI = 0; phenoI < numberRealPheno; ++phenoI) {

			if (variantCorrelationsPruned.rows() > 1) {

				timeStart = System.currentTimeMillis();

				final double geneChi2Sum = geneVariantPhenotypeMatrix.getCol(phenoI).aggregate(DoubleFunctions.plus, DoubleFunctions.square);

				timeStop = System.currentTimeMillis();
				timeInCalculatingRealSumChi2 += (timeStop - timeStart);

				timeStart = System.currentTimeMillis();

				int x = 0;
				int perm = 0;
				int currentNumberPermutations = 0;

				do {

					//LOGGER.debug("Start do with current permutations " + currentNumberPermutations + " x = " + x + " end: " + (currentNumberPermutations + PERMUTATION_STEP) + "?" + ((currentNumberPermutations + PERMUTATION_STEP) < maxNrPermutations));

					if (currentNumberPermutations >= currentNumberPermutationsCalculated) {
						runPermutationsUsingEigenValues(geneChi2SumNull, eigenValues, randomChi2, currentNumberPermutationsCalculated, currentNumberPermutationsCalculated + PERMUTATION_STEP, rnd);
						currentNumberPermutationsCalculated += PERMUTATION_STEP;
					}
					currentNumberPermutations += PERMUTATION_STEP;

					while (++perm < currentNumberPermutations) {
						if (geneChi2Sum < geneChi2SumNull[perm]) {
							x++;
						}
					}

					//LOGGER.debug("test");
				} while (x < 10 && currentNumberPermutations < maxNrPermutations);

				timeStop = System.currentTimeMillis();
				timeInComparingRealChi2ToPermutationChi2 += (timeStop - timeStart);

				timeStart = System.currentTimeMillis();

				double p = (x + 0.5) / (double) (currentNumberPermutations + 1);

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
				genePvalues.setElementQuick(geneI, phenoI, Double.NaN);
				countNoVariants++;
			}

		}

		// Do exactly the same thing but now for the null GWAS
		for (int nullPhenoI = 0; nullPhenoI < NUMBER_RANDOM_PHENO; ++nullPhenoI) {

			if (variantCorrelationsPruned.rows() > 1) {

				timeStart = System.currentTimeMillis();

				final double geneChi2Sum = nullGwasZscores.getCol(nullPhenoI).aggregate(DoubleFunctions.plus, DoubleFunctions.square);

				timeStop = System.currentTimeMillis();
				timeInCalculatingRealSumChi2 += (timeStop - timeStart);

				timeStart = System.currentTimeMillis();

				int x = 0;
				int perm = 0;
				int currentNumberPermutations = 0;

				do {

					//LOGGER.debug("Start do with current permutations " + currentNumberPermutations + " x = " + x + " end: " + (currentNumberPermutations + PERMUTATION_STEP) + "?" + ((currentNumberPermutations + PERMUTATION_STEP) < maxNrPermutations));

					if (currentNumberPermutations >= currentNumberPermutationsCalculated) {
						runPermutationsUsingEigenValues(geneChi2SumNull, eigenValues, randomChi2, currentNumberPermutationsCalculated, currentNumberPermutationsCalculated + PERMUTATION_STEP, rnd);
						currentNumberPermutationsCalculated += PERMUTATION_STEP;
					}
					currentNumberPermutations += PERMUTATION_STEP;

					while (++perm < currentNumberPermutations) {
						if (geneChi2Sum < geneChi2SumNull[perm]) {
							x++;
						}
					}

					//LOGGER.debug("test");
				} while (x < 10 && currentNumberPermutations < maxNrPermutations);

				double p = (x + 0.5) / (double) (currentNumberPermutations + 1);

				if (p == 1) {
					p = 0.99999d;
				}
				if (p < 1e-300d) {
					p = 1e-300d;
				}

				genePvaluesNullGwas.setElementQuick(geneI, nullPhenoI, p);

				timeStop = System.currentTimeMillis();
				timeInCalculatingPvalue += (timeStop - timeStart);

			} else if (variantCorrelationsPruned.rows() == 1) {

				//Always row 0
				double p = ZScores.zToP(-Math.abs(nullGwasZscores.getElementQuick(0, nullPhenoI)));
				if (p == 1) {
					p = 0.99999d;
				}
				if (p < 1e-300d) {
					p = 1e-300d;
				}

				genePvaluesNullGwas.setElementQuick(geneI, nullPhenoI, p);

			} else {
				//no variants in or near gene
				//genePValueDistribution[(int) (20d * 0.99999d)]++;
				genePvaluesNullGwas.setElementQuick(geneI, nullPhenoI, Double.NaN);
			}

		}
		
		geneMaxPermutationCount.setElementQuick(geneI, 0, currentNumberPermutationsCalculated);

	}

	/**
	 * The dosages for each variants will be scaled to have mean of 0 and sd of
	 * 1. This will allow fast correlation calculations
	 *
	 * @param chr
	 * @param start
	 * @param stop
	 * @return
	 */
	private DoubleMatrixDataset<String, String> loadVariantScaledDosageMatrix(String chr, int start, int stop) {

		start = start < 0 ? 0 : start;

		LOGGER.debug("Query genotype data: " + chr + ":" + start + "-" + stop);

		//ArrayList<GeneticVariant> variants = new ArrayList<>(64);
		ArrayList<float[]> variantsDosages = new ArrayList<>(64);
		LinkedHashMap<String, Integer> variantHash = new LinkedHashMap<>(64);

		long timeStart = System.currentTimeMillis();

		int v = 0;
		synchronized (this) {
			for (GeneticVariant variant : referenceGenotypes.getVariantsByRange(chr, start, stop)) {
				//variants.add(variant);
				variantsDosages.add(variant.getSampleDosages());
				variantHash.put(variant.getPrimaryVariantId(), v++);
			}
		}

		DoubleMatrixDataset<String, String> dosageDataset = new DoubleMatrixDataset(sampleHash, variantHash);

		DoubleMatrix2D dosageMatrix = dosageDataset.getMatrix();

		v = 0;
		for (float[] variantDosages : variantsDosages) {
			for (int s = 0; s < variantDosages.length; ++s) {
				dosageMatrix.setQuick(s, v, variantDosages[s]);
			}
			v++;
		}

		long timeStop = System.currentTimeMillis();
		timeInLoadingGenotypeDosages += (timeStop - timeStart);

		LOGGER.debug(" * Variants found in region: " + variantsDosages.size());

		timeStart = System.currentTimeMillis();

		//Inplace normalize rows to mean of 0 and sd of 1 too 
		dosageDataset.normalizeColumns();

		timeStop = System.currentTimeMillis();
		timeInNormalizingGenotypeDosages += (timeStop - timeStart);

		return dosageDataset;

	}

	public static LinkedHashMap<String, Integer> createGeneHashRows(final List<Gene> genes) {
		LinkedHashMap<String, Integer> geneHashRows = new LinkedHashMap<>(genes.size());
		for (int geneI = 0; geneI < genes.size(); ++geneI) {
			geneHashRows.put(genes.get(geneI).getGene(), geneI);
		}
		return geneHashRows;
	}

	public static LinkedHashMap<String, Integer> createHashColsFromList(final List<String> list) {
		LinkedHashMap<String, Integer> hashCols = new LinkedHashMap<>(list.size());
		for (int i = 0; i < list.size(); ++i) {
			hashCols.put(list.get(i), i);
		}
		return hashCols;
	}

	/**
	 * Fills a subset of the null
	 *
	 * @param geneChi2SumNull vector to fill will null chi2 values
	 * @param eigenValues
	 * @param randomChi2 reference distribution
	 * @param start zero based count
	 * @param stop zero based count, so must be < nrPermutation @return
	 */
	public static void runPermutationsUsingEigenValues(final double[] geneChi2SumNull, final double[] eigenValues, double[] randomChi2, final int start, final int stop, final ThreadLocalRandom rnd) {

		final int randomChi2Size = randomChi2.length;

		//LOGGER.debug("start: " + start + " stop: " + stop);

		int i = start;
		for (int perm = start; perm < stop; perm++) {
			double weightedChi2Perm = 0;
			for (int g = 0; g < eigenValues.length; g++) {

				if (i < randomChi2Size) {
					weightedChi2Perm += eigenValues[g] * randomChi2[i++];
				} else {
					double randomZ = rnd.nextGaussian();
					weightedChi2Perm += eigenValues[g] * (randomZ * randomZ);
				}

			}
			geneChi2SumNull[perm] = weightedChi2Perm;
		}

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

	private static DoubleMatrixDataset<String, String> generateRandomNormalizedPheno(LinkedHashMap<String, Integer> sampleHash) {

		LinkedHashMap<String, Integer> phenoHash = new LinkedHashMap<>();

		for (int i = 0; i < NUMBER_RANDOM_PHENO; ++i) {
			phenoHash.put("RanPheno" + (i + 1), i);
		}

		DoubleMatrixDataset<String, String> phenoData = new DoubleMatrixDataset(sampleHash, phenoHash);

		DoubleMatrix2D phenoMatrix = phenoData.getMatrix();

		final int sampleCount = sampleHash.size();

		IntStream.range(0, NUMBER_RANDOM_PHENO).parallel().forEach(pi -> {

			final ThreadLocalRandom rnd = ThreadLocalRandom.current();
			for (int s = 0; s < sampleCount; ++s) {
				phenoMatrix.setQuick(s, pi, rnd.nextGaussian());
			}

		});

		phenoData.normalizeColumns();

		return phenoData;

	}

	private class calculatorThread implements Runnable {

		private final AtomicInteger counter;

		public calculatorThread(AtomicInteger counter) {
			this.counter = counter;
		}

		@Override
		public void run() {

			int i;
			while ((i = counter.getAndIncrement()) < genes.size()) {
				try {
					final double[] geneChi2SumNull = new double[maxNrPermutations];//This will be recyceld
					runGene(i, geneChi2SumNull);
					pb.step();
				} catch (Exception ex) {
					throw new RuntimeException(ex);
				}
			}

		}

	}

}
