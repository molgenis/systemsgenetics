/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.runners;

import java.io.*;
import java.util.*;
import java.util.zip.GZIPInputStream;

import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import me.tongfei.progressbar.ProgressBar;
import nl.systemsgenetics.downstreamer.Downstreamer;
import nl.systemsgenetics.downstreamer.gene.Gene;
import nl.systemsgenetics.downstreamer.io.IoUtils;
import nl.systemsgenetics.downstreamer.runners.options.OptionsModeCoreg;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;
import umcg.genetica.math.PcaColt;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetFastSubsetLoader;
import umcg.genetica.math.stats.PearsonRToZscoreBinned;
import umcg.genetica.math.stats.ZScores;


/**
 *
 * @author patri
 */
public class CoregulationUtilities {

	private static final Logger LOGGER = LogManager.getLogger(CoregulationUtilities.class);
	
	/**
	 * Take a co-regulation matrix and for each gene determine the number of connections that gene has. Reports at
	 * bonferoni significant z-scores, as well as the sumchisqr for each gene.
	 */
	public static void coregInvestigateNetwork(OptionsModeCoreg options) throws Exception {

		int df = options.getNumberSamplesUsedForCor()-2;
		logInfoMem("Hello there.");

		DoubleMatrixDataset<String, String> zScoreMatrix = DoubleMatrixDataset.loadDoubleBinaryData(options.getGwasZscoreMatrixPath());
		logInfoMem("Loaded " + zScoreMatrix.rows() + " rows x " + zScoreMatrix.columns() + " columns");

		if (options.getGeneInfoFile() != null) {
			LinkedHashMap<String, Gene> genes = IoUtils.readGenesMap(options.getGeneInfoFile());
			LOGGER.info("Loaded " + genes.size() + " genes");

			List<String> genesToKeep = zScoreMatrix.getRowObjects();
			LOGGER.info("Read " + genesToKeep.size() + " genes in correlation matrix");
			genesToKeep.retainAll(genes.keySet());

			LOGGER.info("Retained " + genesToKeep.size() + " genes that overlap with --genes file");
			zScoreMatrix = zScoreMatrix.viewSelection(genesToKeep, genesToKeep);

			if (!genes.keySet().containsAll(zScoreMatrix.getHashRows().keySet())) {
				throw new Exception("Not all genes Co-expression matrix are found in gene mapping file");
			}
		}

		if (!zScoreMatrix.getHashRows().keySet().containsAll(zScoreMatrix.getHashCols().keySet())) {
			throw new Exception("Co-expression matrix is not squared with same row and col names");
		}

		final int genesInMatrix = zScoreMatrix.rows();

		// Calculate z-scores back to pearson R
		PearsonRToZscoreBinned zToR = new PearsonRToZscoreBinned(1000000, 12);

		List<String> colnames = new ArrayList<>();
		colnames.add("z_mean");
		colnames.add("z2_sum");
		colnames.add("z2_mean");
		
		colnames.add("z_nominal");
		colnames.add("z_bonf_sig");
		
		colnames.add("r_mean");
		colnames.add("r2_sum");
		colnames.add("r2_mean");
		colnames.add("r2adj_sum");

		double zNominal = Math.abs(ZScores.pToZTwoTailed(0.05));
		double zBonfSig = Math.abs(ZScores.pToZTwoTailed(0.05/genesInMatrix));

		DoubleMatrixDataset<String, String> output = new DoubleMatrixDataset<>(zScoreMatrix.getRowObjects(), colnames);

		ProgressBar pb = new ProgressBar("Coregulation summary stats", genesInMatrix);
		for (int i = 0; i < genesInMatrix; ++i) {

			double curSumZ = 0;
			double curSumZ2 = 0;
			double curZNominal = 0;
			double curZBonfSig = 0;
			double curSumR = 0;
			double curSumR2 = 0;
			double curSumR2adj = 0;

			for(int j =0 ; j < genesInMatrix; j++) {

				double curZ = zScoreMatrix.getElementQuick(i, j);
				curSumZ = curSumZ + curZ;
				curSumZ2 += (curZ * curZ);

				if (Math.abs(curZ) > zNominal) {
					curZNominal ++;
				}

				if (Math.abs(curZ) > zBonfSig) {
					curZBonfSig ++;
				}

				double curR = zToR.lookupRforZscore(curZ);

				if (Double.isNaN(curR)) {
					LOGGER.debug("point");
				}

				curSumR = curSumR + curR;
				double curR2 = curR * curR;
				curSumR2 = curSumR2 + (curR2);

				double curR2adj = curR2 - ((1-curR2) / df);
				curSumR2adj = curSumR2adj + curR2adj;
			}

			output.setElementQuick(i, 0, curSumZ / genesInMatrix);
			output.setElementQuick(i, 1, curSumZ2);
			output.setElementQuick(i, 2, curSumZ2 / genesInMatrix);
			output.setElementQuick(i, 3, curZNominal);
			output.setElementQuick(i, 4, curZBonfSig);
			output.setElementQuick(i, 5, curSumR / genesInMatrix);
			output.setElementQuick(i, 6, curSumR2);
			output.setElementQuick(i, 7, curSumR2 / genesInMatrix);
			output.setElementQuick(i, 8, curSumR2adj);

			pb.step();
		}

		pb.close();
		LOGGER.info("Done, saving output");
		output.save(options.getOutputBasePath() + ".network_properties.tsv");

	}

	private static void logInfoMem(String msg) {
		Downstreamer.logInfoMem(msg, LOGGER);
	}

	/**
	 * Create a gene gene correlation matrix based on a (eigenvector) matrix.
	 *
	 * @param options
	 * @throws Exception
	 */
	public static void coregCorrelateGenes(OptionsModeCoreg options) throws FileNotFoundException, Exception {

		DoubleMatrixDataset<String, String> expressionMatrix;

		if (options.getGeneInfoFile() != null) {
			final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();

			CSVReader reader = null;
			if (options.getGeneInfoFile().getName().endsWith(".gz")) {
				reader = new CSVReaderBuilder((new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(options.getGeneInfoFile())))))).withSkipLines(0).withCSVParser(parser).build();
			} else {
				reader = new CSVReaderBuilder(new BufferedReader(new FileReader(options.getGeneInfoFile()))).withCSVParser(parser).withSkipLines(0).build();
			}

			final HashSet<String> genes = new HashSet<>();
			String[] nextLine;
			while ((nextLine = reader.readNext()) != null) {
				genes.add(nextLine[0]);
			}
			LOGGER.info("Read " + genes.size() + " genes to load");
			if (options.isTrimGeneNames()) {
				LOGGER.info("Note: the genes to load filter is applied before trimming the gene names.");
			}
			if (options.getGwasZscoreMatrixPath().endsWith(".txt") || options.getGwasZscoreMatrixPath().endsWith("txt.gz")) {
				expressionMatrix = DoubleMatrixDataset.loadSubsetOfTextDoubleData(options.getGwasZscoreMatrixPath(), '\t', genes, null);
			} else {

				DoubleMatrixDatasetFastSubsetLoader loader = new DoubleMatrixDatasetFastSubsetLoader(options.getGwasZscoreMatrixPath());
				Set<String> rows = loader.getAllRowIdentifiers();

				genes.retainAll(rows);

				expressionMatrix = loader.loadSubsetOfRowsBinaryDoubleData(genes);
			}

		} else {
			if (options.getGwasZscoreMatrixPath().endsWith(".txt") || options.getGwasZscoreMatrixPath().endsWith("txt.gz")) {
				expressionMatrix = DoubleMatrixDataset.loadDoubleTextData(options.getGwasZscoreMatrixPath(), '\t');
			} else {
				expressionMatrix = DoubleMatrixDataset.loadDoubleBinaryData(options.getGwasZscoreMatrixPath());
			}
		}

		if (options.isTrimGeneNames()) {
			DownstreamerUtilities.trimEnsemblVersionFromRownames(expressionMatrix);
		}

		// Optionally select a subset of columns to use
		String[] cols = options.getColumnsToExtract();
		if (cols != null) {
			Set<String> columnsToExtract = new HashSet<>();
			for (String colname : cols) {
				if (expressionMatrix.getColObjects().contains(colname)) {
					columnsToExtract.add(colname);
				} else {
					LOGGER.warn(colname + " is missing in input matrix, ommiting col in output");
				}
			}

			expressionMatrix = expressionMatrix.viewColSelection(columnsToExtract);
		}

		// Normalize the input data
		if (options.isNormalizeEigenvectors()) {
			expressionMatrix.normalizeRows();
			expressionMatrix.normalizeColumns();
			LOGGER.info("Data row normalized and then column normalized");
		}

		// Calculate the correlation matrix
		LOGGER.info("Loaded expression matrix with " + expressionMatrix.rows() + " genes and " + expressionMatrix.columns() + " observations");
		DoubleMatrixDataset<String, String> corMatrix = expressionMatrix.viewDice().calculateCorrelationMatrix();
		LOGGER.info("Done calculating correlations");

		// Convert Pearson R to Z-scores
		if (options.isConvertRToZscore()) {
			PearsonRToZscoreBinned r2zScore = new PearsonRToZscoreBinned(10000000, expressionMatrix.columns());
			r2zScore.inplaceRToZ(corMatrix);
			LOGGER.info("Converted correlations to Z-scores");
		}

		// Set diagonal of matrix to zero
		for (int i = 0; i < corMatrix.columns(); ++i) {
			corMatrix.setElementQuick(i, i, 0);
		}
		LOGGER.info("Diagnonal set to zero as this might inflate coregulation towards genes in GWAS loci");

		// Save
		LOGGER.info("Saving correlation matrix to: " + options.getOutputBasePath() + ".dat.gz");
		corMatrix.saveBinary(options.getOutputBasePath());
		LOGGER.info("Correlation matrix saved.");

		// Calculate per gene distribution metrics
		LOGGER.info("Calculating per gene distribution metrics");
		DoubleMatrixDataset<String, String> perGeneDistMetrics = DownstreamerUtilities.calculateDistributionMetricsPerRow(corMatrix);
		perGeneDistMetrics.save(options.getOutputBasePath() + ".coregulation.dist.metrics.txt.gz");
		LOGGER.info("Done");

	}

	/**
	 * Converts a matrix of Pearson R values to z-scores. Diagonal of this
	 * matrix is set to zero.
	 *
	 * @param options
	 * @throws FileNotFoundException
	 * @throws Exception
	 */
	public static void coregConvertRtoZscore(OptionsModeCoreg options) throws FileNotFoundException, Exception {
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
	 * Converts a matrix of Pearson R values to z-scores. Diagonal of this
	 * matrix is set to 1.
	 *
	 * @param options
	 * @throws FileNotFoundException
	 * @throws Exception
	 */
	public static void coregConvertZscoreToR(OptionsModeCoreg options) throws FileNotFoundException, Exception {
		DoubleMatrixDataset<String, String> corMatrix;

		if (options.getGwasZscoreMatrixPath().endsWith(".txt") || options.getGwasZscoreMatrixPath().endsWith("txt.gz")) {
			corMatrix = DoubleMatrixDataset.loadDoubleTextData(options.getGwasZscoreMatrixPath(), '\t');
		} else {
			corMatrix = DoubleMatrixDataset.loadDoubleBinaryData(options.getGwasZscoreMatrixPath());
		}

		PearsonRToZscoreBinned tmp = new PearsonRToZscoreBinned((int)20, options.getNumberSamplesUsedForCor());
		tmp.inplaceZToR(corMatrix);
		LOGGER.info("Converted Z-scores to correlations");

		for (int i = 0; i < corMatrix.columns(); ++i) {
			corMatrix.setElementQuick(i, i, 1);
		}

		corMatrix.saveBinary(options.getOutputBasePath());

		LOGGER.info("Correlation matrix saved to: " + options.getOutputBasePath() + ".datg");
	}


	/**
	 * Load a co-regulation matrix and set any gene-gene correlation between
	 * genes closer than 250kb to zero.
	 *
	 * @param options
	 * @throws IOException
	 * @throws Exception
	 */
	public static void coregRemoveLocalGeneCorrelations(OptionsModeCoreg options) throws IOException, Exception {

		LinkedHashMap<String, Gene> genes = IoUtils.readGenesMap(options.getGeneInfoFile());
		LOGGER.info("Loaded " + genes.size() + " genes");

		DoubleMatrixDataset<String, String> corMatrix = DoubleMatrixDataset.loadDoubleBinaryData(options.getGwasZscoreMatrixPath());

		List<String> genesToKeep = corMatrix.getRowObjects();
		LOGGER.info("Read " + genesToKeep.size() + " genes in correlation matrix");
		genesToKeep.retainAll(genes.keySet());
		LOGGER.info("Retained " + genesToKeep.size() + " genes that overlap with --genes file");
		corMatrix = corMatrix.viewSelection(genesToKeep, genesToKeep);

		if (!corMatrix.getHashRows().keySet().containsAll(corMatrix.getHashCols().keySet())) {
			throw new Exception("Co-expression matrix is not squared with same row and col names");
		}

		if (!genes.keySet().containsAll(corMatrix.getHashRows().keySet())) {
			throw new Exception("Not all genes Co-expression matrix are found in gene mapping file");
		}

		final int genesInMatrix = corMatrix.rows();
		final ArrayList<String> geneOrder = corMatrix.getRowObjects();

		int overlappingGenePairs = 0;

		for (int i = 0; i < genesInMatrix; ++i) {

			//diagnoal always 0
			corMatrix.setElementQuick(i, i, 0);

			Gene geneI = genes.get(geneOrder.get(i));

			for (int j = i + 1; j < genesInMatrix; ++j) {

				Gene geneJ = genes.get(geneOrder.get(j));

				if (geneI.withinDistanceOf(geneJ, options.getCisWindow())) {
					corMatrix.setElementQuick(i, j, 0);
					corMatrix.setElementQuick(j, i, 0);
					++overlappingGenePairs;
				}
			}
		}

		LOGGER.info("Identified " + overlappingGenePairs + " overlapping gene-gene pairs within " + options.getCisWindow() + "b.");
		corMatrix.saveBinary(options.getOutputBasePath());

	}

	/**
	 * Run PCA analysis on a binary matrix using PcaColt.
	 *
	 * @param options
	 * @throws IOException
	 */
	public static void coregDoPcaOnBinMatrix(OptionsModeCoreg options) throws IOException {

		final DoubleMatrixDataset<String, String> dataset = DoubleMatrixDataset.loadDoubleBinaryData(options.getGwasZscoreMatrixPath());

		//if debug is enabled keep cov matrix in memory
		PcaColt pcaRes = new PcaColt(dataset, true, true, LOGGER.isDebugEnabled());

		pcaRes.getEigenvectors().save(options.getOutputBasePath() + "_eigenVectors.txt.gz");
		pcaRes.getEigenValues().save(options.getOutputBasePath() + "_eigenValues.txt.gz");
		pcaRes.getPcs().save(options.getOutputBasePath() + "_pcs.txt.gz");
		if (LOGGER.isDebugEnabled()) {
			pcaRes.getCovMatrix().save(options.getOutputBasePath() + "_correlationMatrix.txt.gz");
		}

	}
}
