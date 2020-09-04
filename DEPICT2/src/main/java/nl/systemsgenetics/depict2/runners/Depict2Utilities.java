package nl.systemsgenetics.depict2.runners;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import nl.systemsgenetics.depict2.Depict2;
import nl.systemsgenetics.depict2.Depict2Options;
import nl.systemsgenetics.depict2.Depict2Step2Results;
import nl.systemsgenetics.depict2.Depict2Step3Results;
import nl.systemsgenetics.depict2.gene.Gene;
import nl.systemsgenetics.depict2.io.ExcelWriter;
import nl.systemsgenetics.depict2.io.IoUtils;
import nl.systemsgenetics.depict2.pathway.PathwayDatabase;
import nl.systemsgenetics.depict2.pathway.PathwayEnrichments;
import org.apache.log4j.Logger;
import umcg.genetica.math.PcaColt;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.PearsonRToZscoreBinned;
import umcg.genetica.math.stats.ZScores;

import java.io.*;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetFastSubsetLoader;

/**
 * Collection of runners that handle pre-processing steps for Depict 2 analysis.
 *
 */
public class Depict2Utilities {

	private static final Logger LOGGER = Logger.getLogger(Depict2Utilities.class);

	/**
	 * Create a gene gene correlation matrix based on a (eigenvector) matrix.
	 *
	 * @param options
	 * @throws Exception
	 */
	public static void correlateGenes(Depict2Options options) throws FileNotFoundException, Exception {

		DoubleMatrixDataset<String, String> expressionMatrix;

		if (options.getGeneInfoFile() != null) {
			final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
			final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(options.getGeneInfoFile()))).withCSVParser(parser).withSkipLines(0).build();

			final HashSet<String> genes = new HashSet<>();
			String[] nextLine;
			while ((nextLine = reader.readNext()) != null) {
				genes.add(nextLine[0]);
			}
			LOGGER.info("Read " + genes.size() + " genes to load");
			if (options.getGwasZscoreMatrixPath().endsWith(".txt") || options.getGwasZscoreMatrixPath().endsWith("txt.gz")) {
				expressionMatrix = DoubleMatrixDataset.loadSubsetOfTextDoubleData(options.getGwasZscoreMatrixPath(), '\t', genes, null);
			} else {

				DoubleMatrixDatasetFastSubsetLoader loader = new DoubleMatrixDatasetFastSubsetLoader(options.getGwasZscoreMatrixPath());
				Map<String, Integer> rows = loader.getOriginalRowMap();

				genes.retainAll(rows.keySet());

				expressionMatrix = loader.loadSubsetOfRowsBinaryDoubleData(genes);
			}

		} else {
			if (options.getGwasZscoreMatrixPath().endsWith(".txt") || options.getGwasZscoreMatrixPath().endsWith("txt.gz")) {
				expressionMatrix = DoubleMatrixDataset.loadDoubleTextData(options.getGwasZscoreMatrixPath(), '\t');
			} else {
				expressionMatrix = DoubleMatrixDataset.loadDoubleBinaryData(options.getGwasZscoreMatrixPath());
			}
		}

		if (options.isNormalizeEigenvectors()) {
			expressionMatrix.normalizeRows();
			expressionMatrix.normalizeColumns();
			LOGGER.info("Data row normalized and then column normalized");
		}

		LOGGER.info("Loaded expression matrix with " + expressionMatrix.rows() + " genes and " + expressionMatrix.columns() + " observations");
		DoubleMatrixDataset<String, String> corMatrix = expressionMatrix.viewDice().calculateCorrelationMatrix();
		LOGGER.info("Done calculating correlations");

		if (options.isCorMatrixZscores()) {
			PearsonRToZscoreBinned r2zScore = new PearsonRToZscoreBinned(10000000, expressionMatrix.columns());
			r2zScore.inplaceRToZ(corMatrix);
			LOGGER.info("Converted correlations to Z-scores");
		}

		for (int i = 0; i < corMatrix.columns(); ++i) {
			corMatrix.setElementQuick(i, i, 0);
		}
		LOGGER.info("Diagnonal set to zero as this might inflate coregulation towards genes in GWAS loci");

		corMatrix.saveBinary(options.getOutputBasePath());

		LOGGER.info("Correlation matrix saved to: " + options.getOutputBasePath() + ".dat");

	}

	/**
	 * Converts a matrix of Pearson R values to z-scores. Diagonal of this
	 * matrix is set to zero.
	 *
	 * @param options
	 * @throws FileNotFoundException
	 * @throws Exception
	 */
	public static void convertRtoZscore(Depict2Options options) throws FileNotFoundException, Exception {
		DoubleMatrixDataset<String, String> corMatrix;

		if (options.getGwasZscoreMatrixPath().endsWith(".txt") || options.getGwasZscoreMatrixPath().endsWith("txt.gz")) {
			corMatrix = DoubleMatrixDataset.loadDoubleTextData(options.getGwasZscoreMatrixPath(), '\t');
		} else {
			corMatrix = DoubleMatrixDataset.loadDoubleBinaryData(options.getGwasZscoreMatrixPath());
		}

		PearsonRToZscoreBinned r2zScore = new PearsonRToZscoreBinned(10000000, options.getNumberSamplesUsedForCor());
		r2zScore.inplaceRToZ(corMatrix);
		LOGGER.info("Converted correlations to Z-scores");

		for (int i = 0; i < corMatrix.columns(); ++i) {
			corMatrix.setElementQuick(i, i, 0);
		}
		LOGGER.info("Diagnonal set to zero as this might inflate coregulation towards genes in GWAS loci");

		corMatrix.saveBinary(options.getOutputBasePath());

		LOGGER.info("Correlation matrix saved to: " + options.getOutputBasePath() + ".dat");

	}

	/**
	 * Run PCA analysis on a binary matrix using PcaColt.
	 *
	 * @param options
	 * @throws IOException
	 */
	public static void doPcaOnBinMatrix(Depict2Options options) throws IOException {

		final DoubleMatrixDataset<String, String> dataset = DoubleMatrixDataset.loadDoubleBinaryData(options.getGwasZscoreMatrixPath());

		//if debug is enabled keep cov matrix in memory
		PcaColt pcaRes = new PcaColt(dataset, true, true, LOGGER.isDebugEnabled());

		pcaRes.getEigenvectors().save(options.getOutputBasePath() + "_eigenVectors.txt");
		pcaRes.getEigenValues().save(options.getOutputBasePath() + "_eigenValues.txt");
		pcaRes.getPcs().save(options.getOutputBasePath() + "_pcs.txt");
		if (LOGGER.isDebugEnabled()) {
			pcaRes.getCovMatrix().save(options.getOutputBasePath() + "_correlationMatrix.txt");
		}

	}

	/**
	 * Utility to get the force normalized gene pvalues with tie resolving.
	 * Shares some duplicate code with Depict2MainAnalysis.run2(). This can be
	 * cleaned up in future
	 *
	 * @param options
	 * @throws Exception
	 */
	public static void getNormalizedGwasGenePvalues(Depict2Options options) throws Exception {

		DoubleMatrixDataset<String, String> genePvalues;
		List<Gene> genes;
		DoubleMatrixDataset<String, String> geneVariantCount;
		DoubleMatrixDataset<String, String> geneMaxSnpZscore;

		LOGGER.info("Continuing previous analysis by loading gene p-values");
		if (new File(options.getOutputBasePath() + "_genePvalues.dat").exists()) {
			genePvalues = DoubleMatrixDataset.loadDoubleBinaryData(options.getOutputBasePath() + "_genePvalues");
			// Always load to avoid nullpointers
			geneMaxSnpZscore = DoubleMatrixDataset.loadDoubleBinaryData(options.getOutputBasePath() + "_geneMaxSnpScores");

		} else {
			LOGGER.fatal("Could not find gene pvalues at: " + options.getOutputBasePath() + "_genePvalues.dat");
			LOGGER.fatal("First use --mode RUN to calculate gene p-values");
			return;
		}

		geneVariantCount = DoubleMatrixDataset.loadDoubleTextData(options.getOutputBasePath() + "_geneVariantCount.txt", '\t');
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
		normalizedGwasGeneScores.save(options.getOutputBasePath() + "_normalizedGenePvalues.txt");

		LOGGER.info("Done");

	}

	/**
	 * Re-generate the excel file from existing files in the intermediate
	 * folder.
	 *
	 * @param options
	 * @throws Exception
	 */
	public static void generateExcelFromIntermediates(Depict2Options options) throws Exception {

		Depict2Step2Results step2 = loadExistingStep2Restuls(options);

		ExcelWriter writer = new ExcelWriter(step2.getGenePvalues().getColObjects(), options);

		writer.saveStep2Excel(step2);
		writer.saveGenePvalueExcel(step2.getGenePvalues());

		if (options.getPathwayDatabasesToAnnotateWithGwas().size() >= 1) {
			Depict2Step3Results step3 = Depict2MainAnalysis.step3(options);
			writer.saveStep3Excel(step2, step3);
		}

	}

	public static Depict2Step2Results loadExistingStep2Restuls(Depict2Options options) throws Exception {

		DoubleMatrixDataset<String, String> genePvalues = DoubleMatrixDataset.loadDoubleBinaryData(options.getRun1BasePath() + "_genePvalues");
		final List<PathwayDatabase> pathwayDatabases = options.getPathwayDatabases();
		ArrayList<PathwayEnrichments> pathwayEnrichments = new ArrayList<>(pathwayDatabases.size());
		for (PathwayDatabase pathwayDatabase : pathwayDatabases) {
			pathwayEnrichments.add(new PathwayEnrichments(pathwayDatabase, options.getIntermediateFolder(), options.isExcludeHla()));
		}
		return new Depict2Step2Results(pathwayEnrichments, genePvalues);

	}

}
