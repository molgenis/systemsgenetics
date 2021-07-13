package nl.systemsgenetics.downstreamer.runners;

import cern.colt.function.tdouble.DoubleFunction;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import com.opencsv.*;
import me.tongfei.progressbar.ProgressBar;
import me.tongfei.progressbar.ProgressBarStyle;
import nl.systemsgenetics.downstreamer.DownstreamerOptions;
import nl.systemsgenetics.downstreamer.io.GwasSummStats;
import nl.systemsgenetics.downstreamer.io.IoUtils;
import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetFastSubsetLoader;
import umcg.genetica.math.stats.ZScores;

import java.io.*;
import java.util.*;
import java.util.stream.IntStream;
import java.util.zip.GZIPInputStream;

/**
 * Collection of runners that handle conversion of Depict 2 related files.
 *
 */
public class DownstreamerConverters {

	private static final Logger LOGGER = Logger.getLogger(DownstreamerConverters.class);

	public static void convertTxtToBin(DownstreamerOptions options) throws IOException, Exception {

		if (options.getConversionColumnIncludeFilter() != null && !options.getConversionColumnIncludeFilter().exists()) {
			throw new FileNotFoundException(options.getConversionColumnIncludeFilter().getAbsolutePath() + " (The system cannot find the file specified)");
		}

		if (options.getConversionRowIncludeFilter() != null && !options.getConversionRowIncludeFilter().exists()) {
			throw new FileNotFoundException(options.getConversionRowIncludeFilter().getAbsolutePath() + " (The system cannot find the file specified)");
		}

		final List<String> variantsInZscoreMatrix = DoubleMatrixDataset.readDoubleTextDataRowNames(options.getGwasZscoreMatrixPath(), '\t');
		final List<String> phenotypesInZscoreMatrix = DoubleMatrixDataset.readDoubleTextDataColNames(options.getGwasZscoreMatrixPath(), '\t');

		HashSet<String> phenotypesHashSet = new HashSet<>(phenotypesInZscoreMatrix.size());

		for (String pheno : phenotypesInZscoreMatrix) {
			if (!phenotypesHashSet.add(pheno)) {
				throw new Exception("GWAS matrix contains a duplicate phenotype column: " + pheno);
			}
		}

		HashSet<String> variantsHashSet = new HashSet<>(variantsInZscoreMatrix.size());
		HashSet<String> variantsWithDuplicates = new HashSet<>();

		for (String variant : variantsInZscoreMatrix) {
			if (!variantsHashSet.add(variant)) {
				variantsWithDuplicates.add(variant);
			}
		}

		DoubleMatrixDataset<String, String> matrix;

		if (variantsWithDuplicates.size() > 0) {

			File excludedVariantsFile = new File(options.getOutputBasePath() + "_excludedVariantsDuplicates.txt");

			final CSVWriter excludedVariantWriter = new CSVWriter(new FileWriter(excludedVariantsFile), '\t', '\0', '\0', "\n");
			final String[] outputLine = new String[1];
			outputLine[0] = "ExcludedVariants";
			excludedVariantWriter.writeNext(outputLine);

			for (String dupVariant : variantsWithDuplicates) {
				outputLine[0] = dupVariant;
				excludedVariantWriter.writeNext(outputLine);
			}
			excludedVariantWriter.close();

			LOGGER.info("Found " + variantsWithDuplicates.size() + " duplicate variants, these are excluded from the conversion. For a full list of excluded variants, see: " + excludedVariantsFile.getPath());

			variantsHashSet.removeAll(variantsWithDuplicates);
			matrix = DoubleMatrixDataset.loadSubsetOfTextDoubleData(options.getGwasZscoreMatrixPath(), '\t', variantsHashSet, null);

		} else {
			matrix = DoubleMatrixDataset.loadDoubleTextData(options.getGwasZscoreMatrixPath(), '\t');
		}

		ArrayList<String> allVariants = matrix.getRowObjects();
		ArrayList<String> variantsToExclude = new ArrayList<>();

		if (options.isPvalueToZscore()) {
			DoubleMatrix2D matrixContent = matrix.getMatrix();

			int rows = matrixContent.rows();
			int cols = matrixContent.columns();

			rows:
			for (int r = 0; r < rows; ++r) {
				for (int c = 0; c < cols; ++c) {

					double pvalue = matrixContent.getQuick(r, c);

					if (Double.isNaN(pvalue) || pvalue < 0 || pvalue > 1d) {
						variantsToExclude.add(allVariants.get(c));
						continue rows;
					}

					matrixContent.setQuick(r, c, ZScores.pToZTwoTailed(pvalue));

				}
			}

		} else {
			DoubleMatrix2D matrixContent = matrix.getMatrix();

			int rows = matrixContent.rows();
			int cols = matrixContent.columns();

			rows:
			for (int r = 0; r < rows; ++r) {
				for (int c = 0; c < cols; ++c) {

					double value = matrixContent.getQuick(r, c);

					if (Double.isNaN(value)) {
						variantsToExclude.add(allVariants.get(c));
						continue rows;
					}
				}
			}
		}

		if (variantsToExclude.size() > 0 || options.getConversionRowIncludeFilter() != null) {

			File excludedVariantsFile = new File(options.getOutputBasePath() + "_excludedVariantsNaN.txt");

			final CSVWriter excludedVariantWriter = new CSVWriter(new FileWriter(excludedVariantsFile), '\t', '\0', '\0', "\n");
			final String[] outputLine = new String[1];
			outputLine[0] = "ExcludedVariants";
			excludedVariantWriter.writeNext(outputLine);

			for (String dupVariant : variantsToExclude) {
				outputLine[0] = dupVariant;
				excludedVariantWriter.writeNext(outputLine);
			}
			excludedVariantWriter.close();

			LOGGER.info("Encounterd " + variantsToExclude.size() + " variants with NaN values, these are excluded from the conversion. For a full list of excluded variants, see: " + excludedVariantsFile.getPath());

			HashSet<String> variantsToInclude = new HashSet<>(allVariants);
			variantsToInclude.removeAll(variantsToExclude);

			if (options.getConversionRowIncludeFilter() != null) {
				variantsToInclude.retainAll(new HashSet<>(IoUtils.readMatrixAnnotations(options.getConversionRowIncludeFilter())));
			}

			LOGGER.info("Included variants: " + variantsToInclude.size());

			matrix = matrix.viewRowSelection(variantsToInclude);

		}

		if (options.getConversionColumnIncludeFilter() != null) {
			List<String> colsToSelect = IoUtils.readMatrixAnnotations(options.getConversionColumnIncludeFilter());
			LOGGER.info("Number of selected columns: " + colsToSelect.size());
			matrix = matrix.viewColSelection(colsToSelect);
		}

		matrix.saveBinary(options.getOutputBasePath());

	}

	public static void convertExpressionMatrixToBin(DownstreamerOptions options) throws IOException, Exception {

		if (options.getConversionColumnIncludeFilter() != null && !options.getConversionColumnIncludeFilter().exists()) {
			throw new FileNotFoundException(options.getConversionColumnIncludeFilter().getAbsolutePath() + " (The system cannot find the file specified)");
		}

		if (options.getConversionRowIncludeFilter() != null && !options.getConversionRowIncludeFilter().exists()) {
			throw new FileNotFoundException(options.getConversionRowIncludeFilter().getAbsolutePath() + " (The system cannot find the file specified)");
		}

		DoubleMatrixDataset<String, String> matrix = DoubleMatrixDataset.loadDoubleTextData(options.getGwasZscoreMatrixPath(), '\t');

		if (options.isTrimGeneNames()) {
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

		LOGGER.info("Loaded expression matrix with " + matrix.columns() + " samples and " + matrix.rows() + " genes");

		if (options.getConversionColumnIncludeFilter() != null) {
			List<String> colsToSelect = IoUtils.readMatrixAnnotations(options.getConversionColumnIncludeFilter());
			LOGGER.info("Number of selected columns: " + colsToSelect.size());
			matrix = matrix.viewColSelection(colsToSelect);
		}

		if (options.getConversionRowIncludeFilter() != null) {
			List<String> rowsToSelect = IoUtils.readMatrixAnnotations(options.getConversionRowIncludeFilter());

			LinkedHashMap<String, Integer> allGenes = matrix.getHashRows();
			rowsToSelect.retainAll(allGenes.keySet());

			LOGGER.info("Number of selected rows that exist in matrix: " + rowsToSelect.size());
			matrix = matrix.viewRowSelection(rowsToSelect);
		}

//		matrix.normalizeRows();
//		LOGGER.info("Normalized genes to have mean 0 and sd 1");
		matrix = matrix.createColumnForceNormalDuplicate();
		LOGGER.info("Created force normal of each sample");

		ArrayList<String> rowNames = matrix.getRowObjects();
		ArrayList<String> nonNanRowNames = new ArrayList<>(matrix.rows());

		rows:
		for (int r = 0; r < matrix.rows(); ++r) {
			boolean nonZeroValue = false;
			double e;
			for (int c = 0; c < matrix.columns(); ++c) {
				e = matrix.getElementQuick(r, c);
				if (Double.isNaN(e)) {
					continue rows;
				} else if (e != 0) {
					nonZeroValue = true;
				}
			}
			//if here not NaN values;
			if (nonZeroValue) {
				nonNanRowNames.add(rowNames.get(r));
			}

		}

		if (nonNanRowNames.size() < rowNames.size()) {
			matrix = matrix.viewRowSelection(nonNanRowNames);
			LOGGER.info("Removing " + (rowNames.size() - nonNanRowNames.size()) + " rows with only zero's or NaN after normalizing");
		}

		matrix.saveBinary(options.getOutputBasePath());

	}

	public static void convertBinToTxt(DownstreamerOptions options) throws IOException, Exception {

		DoubleMatrixDataset<String, String> matrix = DoubleMatrixDataset.loadDoubleBinaryData(options.getGwasZscoreMatrixPath());

		String[] cols = options.getColumnsToExtract();

		if (cols != null) {
			Set<String> columnsToExtract = new HashSet<>();

			for (String colname : cols) {
				if (matrix.getColObjects().contains(colname)) {
					columnsToExtract.add(colname);
				} else {
					LOGGER.warn(colname + " is missing in input matrix, ommiting col in output");
				}
			}

			matrix = matrix.viewColSelection(columnsToExtract);
		}

		matrix.save(options.getOutputBasePath() + ".txt");

	}

	public static void convertEqtlToBin(DownstreamerOptions options) throws IOException {

		DoubleMatrixDataset<String, String> matrix = DoubleMatrixDataset.loadTransEqtlExpressionMatrix(options.getGwasZscoreMatrixPath());
		matrix.saveBinary(options.getOutputBasePath());

	}

	/**
	 * Quick util that converts a txt matrix of pvalues to zscores using
	 * ZScores.pToZTwoTailed(pvalue)
	 *
	 * @param options
	 * @throws Exception
	 */
	public static void convertPvalueToZscore(DownstreamerOptions options) throws Exception {

		DoubleMatrixDataset<String, String> matrix = DoubleMatrixDataset.loadDoubleTextData(options.getGwasZscoreMatrixPath(), '\t');

		DoubleMatrix2D matrixContent = matrix.getMatrix();
		LOGGER.info("Read data, converting matrix to zscores");

		// Inplace convert to zscores
		IntStream.range(0, matrix.rows()).parallel().forEach(r -> {
			for (int c = 0; c < matrixContent.columns(); ++c) {
				matrixContent.setQuick(r, c, -ZScores.pToZTwoTailed(matrixContent.getQuick(r, c)));
			}
		});
		LOGGER.info("Done converting matrix to zscores, writing results");

		// Write output
		matrix.save(options.getOutputBasePath() + ".txt");

		LOGGER.info("Done");
	}

	public static void convertGct(String gctFile, String outputMatrixPath) throws Exception {

		final HashSet<String> colToExclude = new HashSet<>();
		colToExclude.add("Description");

		DoubleMatrixDataset<String, String> gtexMedianExp = DoubleMatrixDataset.loadDoubleTextDoubleDataExlcudeCols(gctFile, '\t', colToExclude, 2);

		LOGGER.info("Loaded gtext data of " + gtexMedianExp.rows() + " genes in " + gtexMedianExp.columns() + " tissues");

		HashSet<String> includeGtextGenes = new HashSet<>(gtexMedianExp.getHashRows().keySet());
		HashMap<String, String> ensgToGtexMapping = new HashMap<>();

		for (String gtexGene : gtexMedianExp.getHashRows().keySet()) {

			int indexOfPoint = gtexGene.indexOf('.');
			String gene;
			if (indexOfPoint > 0) {
				gene = gtexGene.substring(0, indexOfPoint);
			} else {
				gene = gtexGene;
			}

			if (ensgToGtexMapping.containsKey(gene)) {
				includeGtextGenes.remove(gtexGene);
				includeGtextGenes.remove(ensgToGtexMapping.get(gene));
				LOGGER.debug("Excluding: " + gene);
			}
			ensgToGtexMapping.put(gene, gtexGene);

		}

		if (includeGtextGenes.size() < gtexMedianExp.rows()) {
			LOGGER.info("Excluding " + (gtexMedianExp.rows() - includeGtextGenes.size()) + " duplicated genes");
			gtexMedianExp = gtexMedianExp.viewRowSelection(includeGtextGenes);
		}

		LinkedHashMap<String, Integer> newHashRow = new LinkedHashMap<>(gtexMedianExp.rows());

		for (Map.Entry<String, Integer> rowEntry : gtexMedianExp.getHashRows().entrySet()) {

			String gene = rowEntry.getKey();

			int indexOfPoint = gene.indexOf('.');
			if (indexOfPoint > 0) {
				gene = gene.substring(0, indexOfPoint);
			}

			newHashRow.put(gene, rowEntry.getValue());

		}

		gtexMedianExp.setHashRows(newHashRow);

		DoubleMatrix2D gtexMedianExpMatrix = gtexMedianExp.getMatrix();
		gtexMedianExpMatrix.assign(new DoubleFunction() {
			@Override
			public double apply(double argument) {
				return Math.log10(argument + 1);
			}
		});

		gtexMedianExp.normalizeRows();

		ArrayList<String> rowNames = gtexMedianExp.getRowObjects();
		ArrayList<String> nonNanRowNames = new ArrayList<>(gtexMedianExp.rows());

		rows:
		for (int r = 0; r < gtexMedianExp.rows(); ++r) {
			for (int c = 0; c < gtexMedianExp.columns(); ++c) {
				if (Double.isNaN(gtexMedianExp.getElementQuick(r, c))) {
					LOGGER.debug("Excluding NaN: " + rowNames.get(r));
					continue rows;
				}
			}
			nonNanRowNames.add(rowNames.get(r));
		}

		if (nonNanRowNames.size() < rowNames.size()) {
			gtexMedianExp = gtexMedianExp.viewRowSelection(nonNanRowNames);
			LOGGER.info("Removing " + (rowNames.size() - nonNanRowNames.size()) + " rows with NaN");
		}

		gtexMedianExp.saveBinary(outputMatrixPath);

	}

	@SuppressWarnings("empty-statement")
	public static void mergeBinMatrix(DownstreamerOptions options) throws IOException, Exception {

		BufferedReader inputReader = new BufferedReader(new FileReader(options.getGwasZscoreMatrixPath()));
		ArrayList<DoubleMatrixDatasetFastSubsetLoader> binMatrices = new ArrayList<>();
		String line;
		while ((line = inputReader.readLine()) != null) {
			if (line.endsWith(".dat")) {
				line = line.substring(0, line.length() - 4);
			}
			binMatrices.add(new DoubleMatrixDatasetFastSubsetLoader(line));
		}

		LinkedHashSet<String> mergedColNames = new LinkedHashSet<>(binMatrices.size());
		LinkedHashSet<String> rowNameIntersection = new LinkedHashSet<>();

		int duplicateColNamesDetected = 0;

		for (DoubleMatrixDatasetFastSubsetLoader datasetLoader : binMatrices) {
			// Put the variant set in memory to avoid having to loop it later on
			if (rowNameIntersection.isEmpty()) {
				rowNameIntersection.addAll(datasetLoader.getOriginalRowMap().keySet());
			} else {
				rowNameIntersection.retainAll(datasetLoader.getOriginalRowMap().keySet());
			}

			for (String newCol : datasetLoader.getOriginalColMap().keySet()) {

				if (mergedColNames.contains(newCol)) {
					int i = 1;
					duplicateColNamesDetected++;
					//intentional empty while loop.
					while (mergedColNames.contains(newCol + "_" + ++i));
					newCol = newCol + "_" + i;
				}
				mergedColNames.add(newCol);
			}

		}

		LOGGER.info("Duplicated col names: " + duplicateColNamesDetected + ". They are included in output file but with a suffix");

		final DoubleMatrixDataset<String, String> mergedData = new DoubleMatrixDataset<>(rowNameIntersection, mergedColNames);

		int mergedCol = 0;
		for (DoubleMatrixDatasetFastSubsetLoader datasetLoader : binMatrices) {
			final DoubleMatrixDataset<String, String> dataset = datasetLoader.loadSubsetOfRowsBinaryDoubleData(rowNameIntersection);

			for (int c = 0; c < dataset.columns(); ++c) {

				mergedData.getCol(mergedCol++).assign(dataset.getCol(c));

			}

		}

		LOGGER.info("Merged data contains: " + mergedData.rows() + " rows and " + mergedData.columns() + " columns");

		mergedData.saveBinary(options.getOutputBasePath());

	}

	public static void mergeConvertTxt(DownstreamerOptions options) throws IOException, Exception {

		if (options.getConversionRowIncludeFilter() != null && !options.getConversionRowIncludeFilter().exists()) {
			throw new FileNotFoundException(options.getConversionRowIncludeFilter().getAbsolutePath() + " (The system cannot find the file specified)");
		}

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(options.getGwasZscoreMatrixPath()))).withCSVParser(parser).build();

		ArrayList<GwasSummStats> gwasSummStats = new ArrayList<>();

		String[] nextLine = reader.readNext();

		if (!"trait".equals(nextLine[0]) || !"file".equals(nextLine[1]) || !"pvalueColumn".equals(nextLine[2])) {
			throw new Exception("Header of file with GWAS summary statistics to use must be: trait<tab>file<tab>pvalueColumn");
		}

		HashSet<String> traits = new HashSet<>();
		while ((nextLine = reader.readNext()) != null) {
			gwasSummStats.add(new GwasSummStats(nextLine[0], nextLine[1], nextLine[2]));
			if (!traits.add(nextLine[0])) {
				throw new Exception("Duplicate trait name: " + nextLine[0]);
			}
		}

		Set<String> overlappingVariants = new HashSet<>();

		try (ProgressBar pb = new ProgressBar("Determining overlapping variants", gwasSummStats.size(), ProgressBarStyle.ASCII)) {
			gwasSummStats.parallelStream().forEach(summStat -> {

				final CSVParser parser2 = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
				try {

					File summStatFile = summStat.getSummStatsFile();

					InputStreamReader isr;
					if (summStatFile.getName().endsWith("gz") || summStatFile.getName().endsWith("bgz")) {
						isr = new InputStreamReader(new GZIPInputStream(new FileInputStream(summStatFile)));
					} else {
						isr = new InputStreamReader(new FileInputStream(summStatFile));
					}

					final CSVReader summStatReader = new CSVReaderBuilder(isr).withSkipLines(1).withCSVParser(parser).build();

					HashSet<String> variantsHashSet = new HashSet<>();
					HashSet<String> variantsWithDuplicates = new HashSet<>();

					String[] summStatNextLine;
					while ((summStatNextLine = summStatReader.readNext()) != null) {
						String variant = summStatNextLine[0];
						if (!variantsHashSet.add(variant)) {
							variantsWithDuplicates.add(variant);
						}
					}

					if (variantsWithDuplicates.size() > 0) {
						File excludedVariantsFile = new File(options.getOutputBasePath() + "_" + summStat.getTrait() + "_excludedVariantsDuplicates.txt");
						final CSVWriter excludedVariantWriter = new CSVWriter(new FileWriter(excludedVariantsFile), '\t', '\0', '\0', "\n");
						final String[] outputLine = new String[1];
						outputLine[0] = "ExcludedVariants";
						excludedVariantWriter.writeNext(outputLine);
						for (String dupVariant : variantsWithDuplicates) {
							outputLine[0] = dupVariant;
							excludedVariantWriter.writeNext(outputLine);
						}
						excludedVariantWriter.close();
						LOGGER.info("Found " + variantsWithDuplicates.size() + " duplicate variants, these are excluded from the conversion. For a full list of excluded variants, see: " + excludedVariantsFile.getPath());

						variantsHashSet.removeAll(variantsWithDuplicates);
					}

					synchronized (overlappingVariants) {
						if (overlappingVariants.isEmpty()) {
							overlappingVariants.addAll(variantsHashSet);
						} else {
							overlappingVariants.retainAll(variantsHashSet);
						}
					}

				} catch (IOException ex) {
					throw new RuntimeException(ex);
				}

				pb.step();

			});
		}

		if (options.getConversionRowIncludeFilter() != null) {
			overlappingVariants.retainAll(new HashSet<>(IoUtils.readMatrixAnnotations(options.getConversionRowIncludeFilter())));
			LOGGER.info("Selected variants that are in the row filter");
		}

		LOGGER.info("Overlapped variants, retained " + overlappingVariants.size());

		ArrayList<String> phenotypes = new ArrayList<>(gwasSummStats.size());
		for (GwasSummStats summStat : gwasSummStats) {
			//not in parralel to keep the order
			phenotypes.add(summStat.getTrait());
		}

		LOGGER.info("Merging over " + phenotypes.size() + " phenotypes");

		// Initialize the output matrix
		DoubleMatrixDataset<String, String> finalMergedPvalueMatrix = new DoubleMatrixDataset<>(overlappingVariants, phenotypes);

		try (ProgressBar pb = new ProgressBar("Meging input data", gwasSummStats.size(), ProgressBarStyle.ASCII)) {

			int i = 0;
			for (GwasSummStats summStat : gwasSummStats) {
				DoubleMatrixDataset<String, String> phenotypeMatrix = summStat.loadSubsetSummStats(overlappingVariants);

				//this is will make sure all variants are in the same order
				phenotypeMatrix = phenotypeMatrix.viewRowSelection(overlappingVariants);

				finalMergedPvalueMatrix.getCol(i++).assign(phenotypeMatrix.getCol(0));
				pb.step();
			}
		}

		ArrayList<String> allVariants = finalMergedPvalueMatrix.getRowObjects();
		ArrayList<String> variantsToExclude = new ArrayList<>();

		if (options.isPvalueToZscore()) {
			DoubleMatrix2D matrixContent = finalMergedPvalueMatrix.getMatrix();

			int rows = matrixContent.rows();
			int cols = matrixContent.columns();

			rows:
			for (int r = 0; r < rows; ++r) {
				for (int c = 0; c < cols; ++c) {

					double pvalue = matrixContent.getQuick(r, c);

					if (Double.isNaN(pvalue) || pvalue < 0 || pvalue > 1d) {
						variantsToExclude.add(allVariants.get(c));
						continue rows;
					}

					matrixContent.setQuick(r, c, ZScores.pToZTwoTailed(pvalue));

				}
			}

		} else {
			DoubleMatrix2D matrixContent = finalMergedPvalueMatrix.getMatrix();

			int rows = matrixContent.rows();
			int cols = matrixContent.columns();

			rows:
			for (int r = 0; r < rows; ++r) {
				for (int c = 0; c < cols; ++c) {

					double value = matrixContent.getQuick(r, c);

					if (Double.isNaN(value)) {
						variantsToExclude.add(allVariants.get(r));
						continue rows;
					}
				}
			}
		}

		if (variantsToExclude.size() > 0) {

			File excludedVariantsFile = new File(options.getOutputBasePath() + "_excludedVariantsNaN.txt");

			final CSVWriter excludedVariantWriter = new CSVWriter(new FileWriter(excludedVariantsFile), '\t', '\0', '\0', "\n");
			final String[] outputLine = new String[1];
			outputLine[0] = "ExcludedVariants";
			excludedVariantWriter.writeNext(outputLine);

			for (String excludedVariant : variantsToExclude) {
				outputLine[0] = excludedVariant;
				excludedVariantWriter.writeNext(outputLine);
			}
			excludedVariantWriter.close();

			LOGGER.info("Encounterd " + variantsToExclude.size() + " variants with NaN values, these are excluded from the conversion. For a full list of excluded variants, see: " + excludedVariantsFile.getPath());

			HashSet<String> variantsToInclude = new HashSet<>(allVariants);
			variantsToInclude.removeAll(variantsToExclude);

			finalMergedPvalueMatrix = finalMergedPvalueMatrix.viewRowSelection(variantsToInclude);

		}

		if (options.getGenotypeBasePath() != null) {
			LOGGER.info("Loading genotype information to convert position  summary statistics to variant IDs");

			File excludedVariantsFile = new File(options.getOutputBasePath() + "_updatedVariantIds.txt");

			final CSVWriter updatedVariantWriter = new CSVWriter(new FileWriter(excludedVariantsFile), '\t', '\0', '\0', "\n");
			final String[] outputLine = new String[2];
			outputLine[0] = "Orginal";
			outputLine[1] = "Updated";
			updatedVariantWriter.writeNext(outputLine);

			RandomAccessGenotypeData genotoypes = IoUtils.loadGenotypes(options, null);

			LinkedHashMap<String, Integer> originalRowHash = finalMergedPvalueMatrix.getHashRows();
			LinkedHashMap<String, Integer> updatedRowHash = new LinkedHashMap<>(originalRowHash.size());

			try (ProgressBar pb = new ProgressBar("Converting variant IDs", originalRowHash.size(), ProgressBarStyle.ASCII)) {
				for (Map.Entry<String, Integer> original : originalRowHash.entrySet()) {

					String originalVariantId = original.getKey();

					String[] splitted = StringUtils.splitPreserveAllTokens(originalVariantId, ':');

					if (splitted.length == 2 || splitted.length == 4) {
						//assuming chr:pos or chr:pos:a1:a2

						String chr = splitted[0];
						int pos = Integer.parseInt(splitted[1]);

						Iterable<GeneticVariant> variantsByPos = genotoypes.getVariantsByPos(chr, pos);

						if (splitted.length == 2) {
							Iterator<GeneticVariant> itt = variantsByPos.iterator();

							if (itt.hasNext()) {
								GeneticVariant variant = itt.next();

								if (itt.hasNext()) {
									//this variant is only variant at position;
									updatedRowHash.put(variant.getPrimaryVariantId(), original.getValue());

									outputLine[0] = originalVariantId;
									outputLine[1] = variant.getPrimaryVariantId();
									updatedVariantWriter.writeNext(outputLine);

								} else {
									//multiple variants, can't match keep original ID
									updatedRowHash.put(original.getKey(), original.getValue());
								}

							} else {
								//No variant at pos
								updatedRowHash.put(original.getKey(), original.getValue());
							}

						} else {
							//splitted.length == 4
							Alleles genotypeGwas = Alleles.createBasedOnString(splitted[2], splitted[3]);

							boolean updated = false;
							variants:
							for (GeneticVariant variant : variantsByPos) {
								if (variant.getVariantAlleles().sameAlleles(genotypeGwas) || (genotypeGwas.isSnp() && variant.getVariantAlleles().sameAlleles(genotypeGwas.getComplement()))) {
									updatedRowHash.put(variant.getPrimaryVariantId(), original.getValue());

									outputLine[0] = originalVariantId;
									outputLine[1] = variant.getPrimaryVariantId();
									updatedVariantWriter.writeNext(outputLine);

									updated = true;
									break variants;
								}
							}
							if (!updated) {
								updatedRowHash.put(original.getKey(), original.getValue());
							}

						}

					} else {
						updatedRowHash.put(original.getKey(), original.getValue());
					}
					pb.step();
				}
			}

			finalMergedPvalueMatrix.setHashRows(updatedRowHash);
			updatedVariantWriter.close();

		}

		finalMergedPvalueMatrix.saveBinary(options.getOutputBasePath());
	}

	public static void tranposeBinMatrix(DownstreamerOptions options) throws IOException {
		DoubleMatrixDataset<String, String> matrix = DoubleMatrixDataset.loadDoubleBinaryData(options.getGwasZscoreMatrixPath());
		matrix = matrix.viewDice();
		matrix.saveBinary(options.getOutputBasePath());
	}

}
