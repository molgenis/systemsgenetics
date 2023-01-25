package nl.systemsgenetics.downstreamer.runners;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.decomposition.DenseDoubleEigenvalueDecomposition;
import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import me.tongfei.progressbar.ProgressBar;
import nl.systemsgenetics.downstreamer.gene.IndexedDouble;
import nl.systemsgenetics.downstreamer.io.BlockDiagonalDoubleMatrixProvider;
import nl.systemsgenetics.downstreamer.io.DoubleMatrixDatasetBlockDiagonalProvider;
import nl.systemsgenetics.downstreamer.runners.options.DownstreamerOptionsDeprecated;
import nl.systemsgenetics.downstreamer.DownstreamerStep2Results;
import nl.systemsgenetics.downstreamer.DownstreamerStep3Results;
import nl.systemsgenetics.downstreamer.gene.Gene;
import nl.systemsgenetics.downstreamer.io.ExcelWriter;
import nl.systemsgenetics.downstreamer.io.IoUtils;
import nl.systemsgenetics.downstreamer.pathway.PathwayEnrichments;
import nl.systemsgenetics.downstreamer.runners.options.OptionsModeCoreg;
import nl.systemsgenetics.downstreamer.runners.options.OptionsTesting;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.inference.MannWhitneyUTest;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;
import umcg.genetica.math.PcaColt;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.PearsonRToZscoreBinned;
import umcg.genetica.math.stats.ZScores;

import java.io.*;
import java.util.*;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.IntStream;
import java.util.zip.GZIPInputStream;

import nl.systemsgenetics.downstreamer.containers.LeadVariant;
import umcg.genetica.collections.ChrPosTreeMap;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetFastSubsetLoader;

/**
 * Collection of runners that handle pre-processing steps for Depict 2 analysis.
 */
public class DownstreamerUtilities {

	private static final Logger LOGGER = LogManager.getLogger(DownstreamerUtilities.class);

	private static HashMap<String, HashMap<String, NearestVariant>> traitGeneDist = null;

	/**
	 * Util to test the block diagonal eigen decomposition. TODO Remove on
	 * release
	 *
	 * @param options
	 * @throws Exception
	 */
	public static void testEigenDecomposition(OptionsTesting options) throws Exception {

		throw new RuntimeException("Compile error");
//		DoubleMatrixDataset<String, String> sigma = DoubleMatrixDataset.loadDoubleData(options.getSigma().getCanonicalPath());
//
//		BlockDiagonalDoubleMatrixProvider provider = new DoubleMatrixDatasetBlockDiagonalProvider(sigma);
//
//		Map<String, List<Integer>> indexMap = new HashMap<>();
//		BufferedReader reader = new BufferedReader(new FileReader(options.getIndex()));
//
//		String line;
//		while ((line = reader.readLine()) != null) {
//			String[] cur = line.split("\t");
//
//			if (indexMap.containsKey(cur[0])) {
//				indexMap.get(cur[0]).add(Integer.parseInt(cur[1]));
//			} else {
//				List<Integer> tmp = new ArrayList<>();
//				tmp.add(Integer.parseInt(cur[1]));
//				indexMap.put(cur[0], tmp);
//			}
//		}
//
//		List<int[]> indices = new ArrayList<>(indexMap.size());
//		for (String key : indexMap.keySet()) {
//			List<Integer> list = indexMap.get(key);
//			int[] array = new int[list.size()];
//			for (int i = 0; i < list.size(); i++) {
//				array[i] = list.get(i);
//			}
//			indices.add(array);
//		}
//		//-----------------------------------------------------------------------------------------------
//
//		List<String> eigenvectorNames = new ArrayList<>(provider.columns());
//		for (int i = 0; i < provider.columns(); i++) {
//			eigenvectorNames.add("V" + i);
//		}
//		List<String> eigenvalueNames = new ArrayList<>(1);
//		eigenvalueNames.add("eigenvalues");
//
//		DoubleMatrixDataset<String, String> U = new DoubleMatrixDataset<>(provider.getRowNames(), eigenvectorNames);
//		DoubleMatrixDataset<String, String> L = new DoubleMatrixDataset<>(eigenvectorNames, eigenvalueNames);
//
//		//-----------------------------------------------------------------------------------------------
//		// Full
//		DenseDoubleEigenvalueDecomposition eigen = new DenseDoubleEigenvalueDecomposition(sigma.getMatrix());
//
//		U.setMatrix(eigen.getV());
//		L.setMatrix(eigen.getRealEigenvalues().reshape(L.rows(), L.columns()));
//
//		U.save(options.getOutputBasePath() + "_ds_eigenvectors.txt");
//		L.save(options.getOutputBasePath() + "_ds_eigenvalues.txt");
//
//		L.setMatrix(eigen.getImagEigenvalues().reshape(L.rows(), L.columns()));
//		L.save(options.getOutputBasePath() + "_ds_eigenvalues_imaginary.txt");
//
//		L.setMatrix(eigen.getRealEigenvalues().reshape(L.rows(), L.columns()));
//
//		//-----------------------------------------------------------------------------------------------
//		// Block diagonal
//		DoubleMatrixDataset<String, String>[] eigenBlock = DownstreamerRegressionEngine.blockDiagonalEigenDecomposition(provider,
//				indices,
//				options.useJblas());
//
//		U = eigenBlock[1];
//		L = eigenBlock[0];
//
//		U.save(options.getOutputBasePath() + "_ds_eigenvectors_blocked.txt");
//		L.save(options.getOutputBasePath() + "_ds_eigenvalues_blocked.txt");
//
//		//-----------------------------------------------------------------------------------------------
	}

	/**
	 * Utility to get the force normalized gene pvalues with tie resolving.
	 * Shares some duplicate code with Depict2MainAnalysis.run2(). This can be
	 * cleaned up in future
	 *
	 * @param options
	 * @throws Exception
	 */
	public static void getNormalizedGwasGenePvalues(DownstreamerOptionsDeprecated options) throws Exception {
		getNormalizedGwasGenePvaluesReturn(options);
	}

	/**
	 * Utility to get the force normalized gene pvalues with tie resolving.
	 * Shares some duplicate code with Depict2MainAnalysis.run2(). This can be
	 * cleaned up in future
	 *
	 * @param options
	 * @throws Exception
	 */
	public static DoubleMatrixDataset<String, String> getNormalizedGwasGenePvaluesReturn(DownstreamerOptionsDeprecated options) throws Exception {

		DoubleMatrixDataset<String, String> genePvalues;
		List<Gene> genes;
		DoubleMatrixDataset<String, String> geneVariantCount;
		DoubleMatrixDataset<String, String> geneMaxSnpZscore;

		LOGGER.info("Continuing previous analysis by loading gene p-values");
		if (new File(options.getRun1BasePath() + "_genePvalues.dat").exists()) {
			genePvalues = DoubleMatrixDataset.loadDoubleBinaryData(options.getRun1BasePath() + "_genePvalues");
			// Always load to avoid nullpointers
			geneMaxSnpZscore = DoubleMatrixDataset.loadDoubleBinaryData(options.getRun1BasePath() + "_geneMaxSnpScores");
		} else {
			LOGGER.fatal("Could not find gene pvalues at: " + options.getRun1BasePath() + "_genePvalues.dat");
			LOGGER.fatal("First use --mode RUN to calculate gene p-values");
			return null;
		}

		String geneVariantCountFile = options.getRun1BasePath() + "_geneVariantCount.txt";
		if (!new File(geneVariantCountFile).canRead()) {
			geneVariantCountFile += ".gz";
			if (!new File(geneVariantCountFile).canRead()) {
				throw new FileNotFoundException("Cannot find file: " + options.getRun1BasePath() + "_geneVariantCount.txt or " + options.getRun1BasePath() + "_geneVariantCount.txt.gz");
			}
		}

		geneVariantCount = DoubleMatrixDataset.loadDoubleTextData(geneVariantCountFile, '\t');
		LOGGER.info("Gene p-values loaded");
		genes = IoUtils.readGenes(options.getGeneInfoFile());
		LOGGER.info("Loaded " + genes.size() + " genes");

		// Identify genes with at least one variant in window
		final HashSet<String> selectedGenes = new HashSet<>();
		final ArrayList<String> allGenes = geneVariantCount.getRowObjects();
		final int totalGeneCount = allGenes.size();
		for (int g = 0; g < totalGeneCount; ++g) {
			if (geneVariantCount.getElementQuick(g, 0) > 0) {
				selectedGenes.add(allGenes.get(g));
			}
		}

		// Select genes that have a gene pvalue, to avoid issues with the normalization, and to keep consistency
		// with the PathwayEnrichments.
		genePvalues = genePvalues.viewRowSelection(selectedGenes);

		LOGGER.info(genePvalues.rows() + " have a gene pvalue");
		final DoubleMatrix2D matrix = genePvalues.getMatrix();

		// Inplace convert gene p-values to z-scores
		IntStream.range(0, matrix.rows()).parallel().forEach(r -> {
			for (int c = 0; c < matrix.columns(); ++c) {
				matrix.setQuick(r, c, -ZScores.pToZTwoTailed(matrix.getQuick(r, c)));
			}
		});

		LOGGER.info("Force normalizing gene p-values / z-scores");
		DoubleMatrixDataset<String, String> normalizedGwasGeneScores;
		normalizedGwasGeneScores = PathwayEnrichments.createColumnForceNormalDuplicate(genePvalues, geneMaxSnpZscore);
		normalizedGwasGeneScores.save(options.getOutputBasePath() + "_normalizedGenePvalues.txt.gz");

		return (normalizedGwasGeneScores);
	}

	/**
	 * Re-generate the excel file from existing files in the intermediate
	 * folder.
	 *
	 * @param options
	 * @throws Exception
	 */
	public static void generateExcelFromIntermediates(DownstreamerOptionsDeprecated options) throws Exception {
throw new RuntimeException("not implemented");
//		DownstreamerStep2Results step2 = IoUtils.loadExistingStep2Results(options);
//		ExcelWriter writer = new ExcelWriter(step2.getGenePvalues().getColObjects(), options);
//
//		writer.saveStep2Excel(step2);
//		//writer.saveGenePvalueExcel(step2.getGenePvalues());
//
//		if (options.getPathwayDatabasesToAnnotateWithGwas().size() >= 1) {
//			DownstreamerStep3Results step3 = DownstreamerMainAnalysis.step3(options);
//			writer.saveStep3Excel(step2, step3);
//		}

	}

	/**
	 * Generate an excel file with the z-scores of the pathways for all bonf.
	 * sig. genes and pathways.
	 *
	 * @param options
	 * @throws Exception
	 */
	public static void generatePathwayLoadingExcel(DownstreamerOptionsDeprecated options) throws Exception {
		throw new RuntimeException("Not implemented");
//		DownstreamerStep2Results step2 = IoUtils.loadExistingStep2Results(options);
//
//		ExcelWriter writer = new ExcelWriter(step2.getGenePvalues().getColObjects(), options);
//		writer.savePathwayLoadings(step2);
	}

	/**
	 * Calculate skewness, kurtosis, mean and variance of the null distribution
	 * for a pathway database
	 *
	 * @throws Exception
	 */
	public static DoubleMatrixDataset<String, String> calculateDistributionMetricsPerRow(DoubleMatrixDataset<String, String> data) {

		List<String> colnames = new ArrayList<>();
		colnames.add("mean");
		colnames.add("sd");
		colnames.add("skewness");
		colnames.add("kurtosis");
		colnames.add("max");
		colnames.add("min");

		DoubleMatrixDataset<String, String> outputStats = new DoubleMatrixDataset<>(data.getRowObjects(), colnames);

		for (int r = 0; r < data.rows(); r++) {
			DescriptiveStatistics curDesc = new DescriptiveStatistics(data.getRow(r).toArray());

			int c = 0;
			outputStats.setElementQuick(r, c++, curDesc.getMean());
			outputStats.setElementQuick(r, c++, curDesc.getStandardDeviation());
			outputStats.setElementQuick(r, c++, curDesc.getSkewness());
			outputStats.setElementQuick(r, c++, curDesc.getKurtosis());
			outputStats.setElementQuick(r, c++, curDesc.getMax());
			outputStats.setElementQuick(r, c++, curDesc.getMin());
		}

		return outputStats;
	}

	/**
	 * Input a normalized expression matrix, and perform a t-test / MannWhitney
	 * between all the samples indicated in the grouping file and the rest. Is
	 * parallelizable.
	 *
	 * @param options
	 */
	public static void generateMarkerGenes(DownstreamerOptionsDeprecated options) throws Exception {

		Set<String> allListedSamples = new HashSet<>();
		Map<String, Set<String>> sampleGroups = new HashMap<>();
		BufferedReader reader = new BufferedReader(new FileReader(options.getX()));
		/*		while (reader.ready()) {
			String[] line = reader.readLine().split("\t");

			Set<String> curSamples = new HashSet<>();
			Collections.addAll(curSamples, Arrays.copyOfRange(line, 1, line.length));

			sampleGroups.put(line[0], curSamples);
			allListedSamples.addAll(curSamples);
		}*/

		String idfile = options.getGwasZscoreMatrixPath() + ".cols.txt";
		if (!new File(idfile).canRead()) {
			idfile += ".gz";
			if (!new File(idfile).canRead()) {
				throw new FileNotFoundException("Cannot find file: " + options.getGwasZscoreMatrixPath() + ".cols.txt or " + options.getGwasZscoreMatrixPath() + ".cols.txt.gz");
			}
		}
		Set<String> availIds = new HashSet<>(DoubleMatrixDataset.readDoubleTextDataRowNames(idfile, '\t'));
		LOGGER.info("Read " + availIds.size() + " available samples");

		while (reader.ready()) {
			String[] line = reader.readLine().split("\t");
			if (line.length != 2) {
				throw new IllegalArgumentException("Line in grouping file contains != 2 columns");
			}

			if (availIds.contains(line[1])) {
				if (sampleGroups.containsKey(line[0])) {
					sampleGroups.get(line[0]).add(line[1]);
				} else {
					Set<String> curSet = new HashSet<>();
					curSet.add(line[1]);
					sampleGroups.put(line[0], curSet);
				}
				allListedSamples.add(line[1]);
			}

		}

		LOGGER.info("Read " + sampleGroups.size() + " sample groups over " + allListedSamples.size() + " samples");

		List<Gene> genes = IoUtils.readGenes(options.getGeneInfoFile());
		Set<String> geneIds = new HashSet<>();
		for (Gene gene : genes) {
			geneIds.add(gene.getGene());
		}

		LOGGER.info("Read " + genes.size() + " genes groups to test");

		// Bulk load expression matrix into memory, could convert to streaming algorithm for more mem efficiency
		LOGGER.info("Loading expression data");

		String rowIdFile = options.getGwasZscoreMatrixPath() + ".rows.txt";
		if (!new File(rowIdFile).canRead()) {
			rowIdFile += ".gz";
			if (!new File(rowIdFile).canRead()) {
				throw new FileNotFoundException("Cannot find file: " + options.getGwasZscoreMatrixPath() + ".rows.txt or " + options.getGwasZscoreMatrixPath() + ".rows.txt.gz");
			}
		}
		geneIds.retainAll(DoubleMatrixDataset.readDoubleTextDataRowNames(rowIdFile, '\t'));

		final DoubleMatrixDataset<String, String> expression = DoubleMatrixDataset.loadSubsetOfRowsBinaryDoubleData(options.getGwasZscoreMatrixPath(), geneIds).viewColSelection(allListedSamples);
		LOGGER.info("Done Loading expression data");

		// Ttest object, and storage for output
		//final TTest tTest = new TTest();
		// Mann Whitney test
		final MannWhitneyUTest mannWhitneyTest = new MannWhitneyUTest();

		final DoubleMatrixDataset<String, String> pvalues = new DoubleMatrixDataset<>(geneIds, sampleGroups.keySet());
		final DoubleMatrixDataset<String, String> stats = new DoubleMatrixDataset<>(geneIds, sampleGroups.keySet());
		final DoubleMatrixDataset<String, String> deltas = new DoubleMatrixDataset<>(geneIds, sampleGroups.keySet());
		final DoubleMatrixDataset<String, String> meansA = new DoubleMatrixDataset<>(geneIds, sampleGroups.keySet());
		final DoubleMatrixDataset<String, String> meansB = new DoubleMatrixDataset<>(geneIds, sampleGroups.keySet());

		//LOGGER.info("Running T-tests");
		LOGGER.info("Running MannWhitney-tests");
		final ProgressBar bp = new ProgressBar("genes", geneIds.size());
		bp.step();
		// Loop over genes, to make more memory efficient, could load one gene at the time
		// Currently with expression matrix of 19x17k its quite doable
		//for (String gene : geneIds) {

		ForkJoinPool forkJoinPool = null;
		try {
			forkJoinPool = new ForkJoinPool(options.getNumberOfThreadsToUse());
			forkJoinPool.submit(() -> geneIds.parallelStream().forEach(gene -> {
				Set<String> curGene = new HashSet<>();
				curGene.add(gene);

				// Loop over the differnet tissue groups
				for (String tissue : sampleGroups.keySet()) {

					// Determine the groups to test
					Set<String> groupA = sampleGroups.get(tissue);
					Set<String> groupB = new HashSet<>();

					for (String tissue2 : sampleGroups.keySet()) {
						if (!tissue2.equals(tissue)) {
							groupB.addAll(sampleGroups.get(tissue2));
						}
					}

					// Subset the 2 sets
					DoubleMatrix1D groupAValues = expression.viewSelection(curGene, groupA).getRow(gene);
					DoubleMatrix1D groupBValues = expression.viewSelection(curGene, groupB).getRow(gene);

					// Run T-test
					//double tstat = tTest.t(groupAValues.toArray(), groupBValues.toArray());
					//double pval = tTest.tTest(groupAValues.toArray(), groupBValues.toArray());
					// Mann Whitney U
					double[] groupAValueArray = groupAValues.toArray();
					double[] groupBValueArray = groupBValues.toArray();

					double pval = mannWhitneyTest.mannWhitneyUTest(groupAValueArray, groupBValueArray);
					double stat = mannWhitneyTest.mannWhitneyU(groupAValueArray, groupBValueArray);
					double meanA = mean(groupAValueArray);
					double meanB = mean(groupBValueArray);

					// Save to output matrix
					pvalues.setElement(gene, tissue, pval);
					stats.setElement(gene, tissue, stat);
					meansA.setElement(gene, tissue, meanA);
					meansB.setElement(gene, tissue, meanB);
					deltas.setElement(gene, tissue, meanA - meanB);

				}

				bp.step();
			})).get();

		} catch (Exception e) {
			throw new RuntimeException(e);
		} finally {
			if (forkJoinPool != null) {
				forkJoinPool.shutdown();
				bp.close();
				//LOGGER.info("Done calculating T-tests");
				LOGGER.info("Done calculating MannWhitney-tests");
			}
		}

		// Save output
		pvalues.save(options.getOutputBasePath() + ".mwu.pvalues.txt.gz");
		stats.save(options.getOutputBasePath() + ".mwu.stats.txt.gz");
		deltas.save(options.getOutputBasePath() + ".mwu.delta.txt.gz");
		meansA.save(options.getOutputBasePath() + ".mwu.meanA.txt.gz");
		meansB.save(options.getOutputBasePath() + ".mwu.meanB.txt.gz");

		LOGGER.info("Done");

	}

	/**
	 * HashMap<Trait, Map<Gene, distance> <0 for other chr or outside cis window
	 *
	 * @param options
	 * @return
	 * @throws IOException
	 */
	public static HashMap<String, HashMap<String, NearestVariant>> getDistanceGeneToTopCisSnpPerTrait(
			final DownstreamerOptionsDeprecated options) throws IOException {

		if (traitGeneDist != null) {
			return traitGeneDist;
		}

		Map<String, ChrPosTreeMap<LeadVariant>> indepVariantsAsSummaryStatisticsRecord = IoUtils.loadLeadVariantsPerTrait(options);
		LinkedHashMap<String, Gene> genes = IoUtils.readGenesMap(options.getGeneInfoFile());
		final int cisExtent = options.getCisWindowExtend();

		HashMap<String, HashMap<String, NearestVariant>> traitGeneDist2 = new HashMap<>(indepVariantsAsSummaryStatisticsRecord.size());

		for (Map.Entry<String, ChrPosTreeMap<LeadVariant>> traitEntry : indepVariantsAsSummaryStatisticsRecord.entrySet()) {

			String trait = traitEntry.getKey();
			ChrPosTreeMap<LeadVariant> topHits = traitEntry.getValue();

			HashMap<String, NearestVariant> geneDist = new HashMap<>(genes.size());
			traitGeneDist2.put(trait, geneDist);

			for (Gene gene : genes.values()) {

				int geneStart = Math.min(gene.getStart(), gene.getEnd());
				int geneEnd = Math.max(gene.getStart(), gene.getEnd());

				int minDist = Integer.MAX_VALUE;
				LeadVariant nearestVariant = null;

				for (LeadVariant cisVariant : topHits.getChrRange(gene.getContig(), geneStart - cisExtent, true, geneEnd + cisExtent, true).values()) {
					if (cisVariant.getPos() >= geneStart && cisVariant.getPos() <= geneEnd) {
						minDist = 0;
						nearestVariant = cisVariant;
						continue;
					}

					int dist;
					if (cisVariant.getPos() < geneStart) {
						dist = geneStart - cisVariant.getPos();
					} else {
						dist = cisVariant.getPos() - geneEnd;
					}
					if (dist < minDist) {
						minDist = dist;
						nearestVariant = cisVariant;
					}

				}

				if (minDist == Integer.MAX_VALUE) {
					//no variant found in window
					minDist = -9;
				}

				geneDist.put(gene.getGene(), new NearestVariant(nearestVariant, minDist));

			}

		}

		traitGeneDist = traitGeneDist2;

		return traitGeneDist2;

	}

	/**
	 * Calculate the benjamini hochberg adjusted p-values form a
	 * doublematrixdataset. Preserves the order of the orginal input. Adapted
	 * from:
	 * https://github.com/cBioPortal/cbioportal/blob/master/core/src/main/java/org/mskcc/cbio/portal/stats/BenjaminiHochbergFDR.java
	 * and
	 * https://stats.stackexchange.com/questions/238458/whats-the-formula-for-the-benjamini-hochberg-adjusted-p-value
	 *
	 * @param pvalues
	 * @return
	 */
	public static DoubleMatrixDataset<String, String> adjustPvaluesBenjaminiHochberg(DoubleMatrixDataset<String, String> rawPvalues) {
		DoubleMatrixDataset<String, String> output = new DoubleMatrixDataset<>(rawPvalues.getRowObjects(), rawPvalues.getColObjects());
		int m = rawPvalues.rows();

		for (int c = 0; c < rawPvalues.columns(); c++) {

			String colname = output.getColObjects().get(c);

			// Sort the p-values preserving the ids
			List<DoubleElement> sortedPvalues = new ArrayList<>(rawPvalues.rows());
			for (int r = 0; r < rawPvalues.rows(); r++) {
				sortedPvalues.add(new DoubleElement(rawPvalues.getElementQuick(r, c), rawPvalues.getRowObjects().get(r)));
			}
			sortedPvalues.sort(Comparator.comparing(DoubleElement::getValue));

			List<DoubleElement> adjustedPvalues = new ArrayList<>(sortedPvalues);

			// iterate through all p-values:  largest to smallest
			for (int i = m - 1; i >= 0; i--) {
				if (i == m - 1) {
					adjustedPvalues.set(i, sortedPvalues.get(i));
				} else {
					double unadjustedPvalue = sortedPvalues.get(i).value;
					int divideByM = i + 1;
					double left = adjustedPvalues.get(i + 1).value;
					double right = (m / (double) divideByM) * unadjustedPvalue;
					adjustedPvalues.set(i, new DoubleElement(Math.min(left, right), sortedPvalues.get(i).id));
				}
			}

			for (DoubleElement curElement : adjustedPvalues) {
				output.setElement(curElement.id, colname, curElement.value);
			}

		}

		return output;
	}

	public static double mean(double[] m) {
		double sum = 0;
		for (int i = 0; i < m.length; i++) {
			sum += m[i];
		}
		return sum / m.length;
	}

	public static class NearestVariant {

		private final LeadVariant nearestVariant;
		private final int distance;

		public NearestVariant(LeadVariant nearestVariant, int distance) {
			this.nearestVariant = nearestVariant;
			this.distance = distance;
		}

		public LeadVariant getNearestVariant() {
			return nearestVariant;
		}

		public int getDistance() {
			return distance;
		}

	}

	/**
	 * Links a double with an ID. Used in adjustPvaluesBenjaminiHochberg.
	 */
	private static class DoubleElement {

		protected double value;
		protected String id;

		public DoubleElement(double value, String id) {
			this.value = value;
			this.id = id;
		}

		public double getValue() {
			return value;
		}

		public String getId() {
			return id;
		}
	}

	/**
	 * Trims the ENSG00000000.X from ensembl gene names. Made this a seperate
	 * function as code was duplicated 3 times
	 *
	 * @param matrix
	 * @throws Exception
	 */
	public static void trimEnsemblVersionFromRownames(DoubleMatrixDataset<String, String> matrix) throws Exception {

		LinkedHashMap<String, Integer> oldHash = matrix.getHashRows();
		LinkedHashMap<String, Integer> newHash = new LinkedHashMap<>(oldHash.size());

		for (Map.Entry<String, Integer> oldEntry : oldHash.entrySet()) {

			String oldGeneName = oldEntry.getKey();
			int indexOfPoint = oldGeneName.indexOf('.');

			if (indexOfPoint < 0) {
				if (newHash.put(oldGeneName, oldEntry.getValue()) != null) {
					throw new Exception("Can't trim gene names if this causes duplicate genes: " + oldGeneName);
				}
			} else {
				if (newHash.put(oldGeneName.substring(0, indexOfPoint), oldEntry.getValue()) != null) {
					throw new Exception("Can't trim gene names if this causes duplicate genes: " + oldGeneName);
				}
			}
		}

		matrix.setHashRows(newHash);
	}

}
