package nl.systemsgenetics.depict2;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import com.opencsv.CSVWriter;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.concurrent.ThreadLocalRandom;
import java.util.stream.IntStream;
import nl.systemsgenetics.depict2.development.ExtractCol;
import nl.systemsgenetics.depict2.development.First1000qtl;
import org.apache.commons.cli.ParseException;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.lang.time.DurationFormatUtils;
import org.apache.log4j.ConsoleAppender;
import org.apache.log4j.FileAppender;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.SimpleLayout;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.multipart.IncompatibleMultiPartGenotypeDataException;
import org.molgenis.genotype.sampleFilter.SampleFilter;
import org.molgenis.genotype.sampleFilter.SampleIdIncludeFilter;
import org.molgenis.genotype.tabix.TabixFileNotFoundException;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variantFilter.VariantCombinedFilter;
import org.molgenis.genotype.variantFilter.VariantFilter;
import org.molgenis.genotype.variantFilter.VariantFilterMaf;
import org.molgenis.genotype.variantFilter.VariantIdIncludeFilter;
import umcg.genetica.math.PcaColt;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetFastSubsetLoader;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.math.stats.PearsonRToZscoreBinned;

/**
 *
 * @author Patrick Deelen
 */
public class Depict2 {

	public static final DecimalFormat LARGE_INT_FORMAT = new DecimalFormat("###,###");
	public static final String VERSION = ResourceBundle.getBundle("verion").getString("application.version");
	private static final DateFormat DATE_TIME_FORMAT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
	private static final Logger LOGGER = Logger.getLogger(Depict2.class);
	private static final String HEADER
			= "  /---------------------------------------\\\n"
			+ "  |                DEPICT2                |\n"
			+ "  |                                       |\n"
			+ "  |  University Medical Center Groningen  |\n"
			+ "  \\---------------------------------------/";

	/**
	 * @param args the command line arguments
	 * @throws java.lang.InterruptedException
	 */
	public static void main(String[] args) throws InterruptedException {

		System.out.println(HEADER);
		System.out.println();
		System.out.println("          --- Version: " + VERSION + " ---");
		System.out.println();
		System.out.println("More information: http://molgenis.org/systemsgenetics");
		System.out.println();

		Date currentDataTime = new Date();
		String startDateTime = DATE_TIME_FORMAT.format(currentDataTime);
		System.out.println("Current date and time: " + startDateTime);
		System.out.println();

		System.out.flush(); //flush to make sure header is before errors
		Thread.sleep(25); //Allows flush to complete

		Depict2Options options;

		if (args.length == 0) {
			Depict2Options.printHelp();
			return;
		}

		try {
			options = new Depict2Options(args);
		} catch (ParseException ex) {
			System.err.println("Error parsing commandline: " + ex.getMessage());
			Depict2Options.printHelp();
			return;
		}

		if (options.getLogFile().getParentFile() != null && !options.getLogFile().getParentFile().isDirectory()) {
			if (!options.getLogFile().getParentFile().mkdirs()) {
				System.err.println("Failed to create output folder: " + options.getLogFile().getParent());
				System.exit(1);
			}
		}

		if (new File(options.getOutputBasePath()).isDirectory()) {
			System.err.println("Specified output path is a directory. Please include a prefix for the output files.");
			return;
		}

		try {
			FileAppender logFileAppender = new FileAppender(new SimpleLayout(), options.getLogFile().getCanonicalPath(), options.getMode() == Depict2Mode.RUN2);
			ConsoleAppender logConsoleInfoAppender = new ConsoleAppender(new InfoOnlyLogLayout());
			Logger.getRootLogger().removeAllAppenders();
			Logger.getRootLogger().addAppender(logFileAppender);

			LOGGER.info("DEPICT" + VERSION);
			LOGGER.info("Current date and time: " + startDateTime);

			Logger.getRootLogger().addAppender(logConsoleInfoAppender);

			if (options.isDebugMode()) {
				Logger.getRootLogger().setLevel(Level.DEBUG);
				options.getDebugFolder().mkdir();
			} else {
				Logger.getRootLogger().setLevel(Level.INFO);
			}

		} catch (IOException e) {
			System.err.println("Failed to create logger: " + e.getMessage());
			System.exit(1);
		}

		options.printOptions();

		try {
			switch (options.getMode()) {
				case CONVERT_EQTL:
					convertEqtlToBin(options);
					break;
				case CONVERT_TXT:
					convertTxtToBin(options);
					break;
				case CONVERT_BIN:
					convertBinToTxt(options);
					break;
				case CONVERT_GTEX:
					ConvertGtexGct.convertGct(options.getGwasZscoreMatrixPath(), options.getOutputBasePath());
					break;
				case CONVERT_EXP:
					convertExpressionMatrixToBin(options);
					break;
				case CONVERT_TXT_MERGE:
					mergeConvertTxt(options);
					break;
				case MERGE_BIN:
					mergeBinMatrix(options);
					break;
				case FIRST1000:
					First1000qtl.printFirst1000(options);
					break;
				case RUN:
					run(options);
					break;
				case RUN2:
					run2(options, null, null, null, null);
					break;
				case CORRELATE_GENES:
					correlateGenes(options);
					break;
				case TRANSPOSE:
					tranposeBinMatrix(options);
					break;
				case PCA:
					doPcaOnBinMatrix(options);
					break;
				case SPECIAL:
					ExtractCol.extract(options.getGwasZscoreMatrixPath(), "GO:0001501", options.getOutputBasePath());
			}
		} catch (TabixFileNotFoundException e) {
			System.err.println("Problem running mode: " + options.getMode());
			System.err.println("Tabix file not found for input data at: " + e.getPath() + "\n"
					+ "Please see README on how to create a tabix file");
			LOGGER.fatal("Tabix file not found for input data at: " + e.getPath(), e);
			System.exit(1);
		} catch (IOException e) {
			System.err.println("Problem running mode: " + options.getMode());
			System.err.println("Error accessing input data: " + e.getMessage());
			System.err.println("See log file for stack trace");
			LOGGER.fatal("Error accessing input data: " + e.getMessage(), e);
			System.exit(1);
		} catch (IncompatibleMultiPartGenotypeDataException e) {
			System.err.println("Problem running mode: " + options.getMode());
			System.err.println("Error combining the impute genotype data files: " + e.getMessage());
			System.err.println("See log file for stack trace");
			LOGGER.fatal("Error combining the impute genotype data files: " + e.getMessage(), e);
			System.exit(1);
		} catch (GenotypeDataException e) {
			System.err.println("Problem running mode: " + options.getMode());
			System.err.println("Error reading genotype data: " + e.getMessage());
			System.err.println("See log file for stack trace");
			LOGGER.fatal("Error reading input data: " + e.getMessage(), e);
			System.exit(1);
		} catch (Exception e) {
			System.err.println("Problem running mode: " + options.getMode());
			System.err.println("Error message: " + e.getMessage());
			System.err.println("See log file for stack trace");
			LOGGER.fatal("Error: " + e.getMessage(), e);
			if (LOGGER.isDebugEnabled()) {
				e.printStackTrace();
			}
			System.exit(1);
		}
		LOGGER.info("Analysis completed");

		currentDataTime = new Date();
		LOGGER.info("Current date and time: " + DATE_TIME_FORMAT.format(currentDataTime));

	}

	public static void run(Depict2Options options) throws IOException, Exception {

		//Test here to prevent chrash after running for a while
		if (!new File(options.getGwasZscoreMatrixPath() + ".dat").exists()) {
			throw new FileNotFoundException("GWAS matrix does not exist at: " + options.getGwasZscoreMatrixPath() + ".dat");
		}

		final List<String> variantsInZscoreMatrix = readMatrixAnnotations(new File(options.getGwasZscoreMatrixPath() + ".rows.txt"));
		final List<String> phenotypesInZscoreMatrix = readMatrixAnnotations(new File(options.getGwasZscoreMatrixPath() + ".cols.txt"));

		LOGGER.info("Number of phenotypes in GWAS matrix: " + LARGE_INT_FORMAT.format(phenotypesInZscoreMatrix.size()));
		LOGGER.info("Number of variants in GWAS matrix: " + LARGE_INT_FORMAT.format(variantsInZscoreMatrix.size()));

		if (options.getVariantFilterFile() != null) {
			HashSet<String> variantsToInclude = readVariantFilterFile(options.getVariantFilterFile());

			Iterator<String> variantsInZscoreMatrixIt = variantsInZscoreMatrix.iterator();
			while (variantsInZscoreMatrixIt.hasNext()) {
				String variant = variantsInZscoreMatrixIt.next();
				if (!variantsToInclude.contains(variant)) {
					variantsInZscoreMatrixIt.remove();
				}

			}
			LOGGER.info("Number of variants after filtering on selected variants: " + LARGE_INT_FORMAT.format(variantsInZscoreMatrix.size()));
		}

		RandomAccessGenotypeData referenceGenotypeData = loadGenotypes(options, variantsInZscoreMatrix);

		LOGGER.info("Done loading genotype data");

		List<Gene> genes = readGenes(options.getGeneInfoFile());

		LOGGER.info("Loaded " + genes.size() + " genes");

		double[] randomChi2 = generateRandomChi2(options.getNumberOfPermutationsRescue(), 500);

		LOGGER.info("Prepared reference null distribution with " + LARGE_INT_FORMAT.format(randomChi2.length) + " values");

		File usedVariantsPerGeneFile = options.isSaveUsedVariantsPerGene() ? new File(options.getOutputBasePath() + "_usedVariantsPerGene.txt") : null;

		GenePvalueCalculator gpc = new GenePvalueCalculator(options.getGwasZscoreMatrixPath(), referenceGenotypeData, genes, options.getWindowExtend(), options.getMaxRBetweenVariants(), options.getNumberOfPermutations(), options.getNumberOfPermutationsRescue(), options.getOutputBasePath(), randomChi2, options.correctForLambdaInflation(), options.getPermutationGeneCorrelations(), options.getPermutationPathwayEnrichment(), options.getDebugFolder(), options.getVariantGeneLinkingFile(), usedVariantsPerGeneFile);

		DoubleMatrixDataset<String, String> genePvalues = gpc.getGenePvalues();
		DoubleMatrixDataset<String, String> genePvaluesNullGwas = gpc.getGenePvaluesNullGwas();
		DoubleMatrixDataset<String, String> geneVariantCount = gpc.getGeneVariantCount();

		genePvalues.saveBinary(options.getOutputBasePath() + "_genePvalues");
		genePvaluesNullGwas.saveBinary(options.getOutputBasePath() + "_genePvaluesNullGwas");
		geneVariantCount.save(options.getOutputBasePath() + "_geneVariantCount.txt");
		if (LOGGER.isDebugEnabled()) {
			gpc.getGeneMaxPermutationCount().save(options.getOutputBasePath() + "_geneMaxPermutationUsed.txt");
			gpc.getGeneRuntime().save(options.getOutputBasePath() + "_geneRuntime.txt");
		}
		LOGGER.info("Gene p-values saved. If needed the analysis can be resummed from this point using --mode RUN2 and exactly the same output path and genes file");

		if (options.getPathwayDatabases().isEmpty()) {
			LOGGER.info("The analysis will now stop since no pathway databases are provided. Use --mode RUN2 and exactly the same output path and genes file to continue");
		} else {
			run2(options, genePvalues, genePvaluesNullGwas, genes, geneVariantCount);
		}
	}

	/**
	 *
	 * @param options
	 * @param genePvalues
	 * @param genePvaluesNullGwas
	 * @param genes
	 * @throws IOException
	 * @throws Exception
	 */
	private static void run2(Depict2Options options, DoubleMatrixDataset<String, String> genePvalues, DoubleMatrixDataset<String, String> genePvaluesNullGwas, List<Gene> genes, DoubleMatrixDataset<String, String> geneVariantCount) throws IOException, Exception {

		options.getIntermediateFolder().mkdir();

		if (options.getMode() == Depict2Mode.RUN2) {
			LOGGER.info("Continuing previous analysis by loading gene p-values");
			if (new File(options.getOutputBasePath() + "_genePvalues.dat").exists()) {
				genePvalues = DoubleMatrixDataset.loadDoubleBinaryData(options.getOutputBasePath() + "_genePvalues");
				genePvaluesNullGwas = DoubleMatrixDataset.loadDoubleBinaryData(options.getOutputBasePath() + "_genePvaluesNullGwas");
			} else if (new File(options.getOutputBasePath() + "_genePvalues.txt").exists()) {
				//This is for some legacy results. New versions of DEPICT2 will not create these files
				genePvalues = DoubleMatrixDataset.loadDoubleTextData(options.getOutputBasePath() + "_genePvalues.txt", '\t');
				genePvaluesNullGwas = DoubleMatrixDataset.loadDoubleTextData(options.getOutputBasePath() + "_genePvaluesNullGwas.txt", '\t');
			} else {
				LOGGER.fatal("Could not find gene pvalues at: " + options.getOutputBasePath() + "_genePvalues.dat");
				LOGGER.fatal("First use --mode RUN to calculate gene p-values");
				return;
			}
			geneVariantCount = DoubleMatrixDataset.loadDoubleTextData(options.getOutputBasePath() + "_geneVariantCount.txt", '\t');
			LOGGER.info("Gene p-values loaded");
			genes = readGenes(options.getGeneInfoFile());
			LOGGER.info("Loaded " + genes.size() + " genes");
		}

		//Identify genes with atleast one variant in window
		final HashSet<String> selectedGenes = new HashSet<>();
		final ArrayList<String> allGenes = geneVariantCount.getRowObjects();
		final int totalGeneCount = allGenes.size();
		for (int g = 0; g < totalGeneCount; ++g) {
			if (geneVariantCount.getElementQuick(g, 0) > 0) {
				selectedGenes.add(allGenes.get(g));
			}
		}

		final DoubleMatrix2D matrix = genePvalues.getMatrix();

		//Inplace convert gene p-values to z-scores
		IntStream.range(0, matrix.rows()).parallel().forEach(r -> {
			for (int c = 0; c < matrix.columns(); ++c) {
				matrix.setQuick(r, c, -ZScores.pToZTwoTailed(matrix.getQuick(r, c)));
			}
		});

		DoubleMatrix2D matrixNull = genePvaluesNullGwas.getMatrix();

		IntStream.range(0, matrixNull.rows()).parallel().forEach(r -> {
			for (int c = 0; c < matrixNull.columns(); ++c) {
				matrixNull.setQuick(r, c, -ZScores.pToZTwoTailed(matrixNull.getQuick(r, c)));
			}
		});

		LOGGER.info("Number of genes with atleast one variant in specified window: " + LARGE_INT_FORMAT.format(selectedGenes.size()));

		final HashSet<String> hlaGenes;
		if (options.isExcludeHla()) {
			hlaGenes = new HashSet<>();
			for (Gene gene : genes) {
				if (gene.getChr().equals("6") && ((gene.getStart() > 20000000 && gene.getStart() < 40000000) || (gene.getStop() > 20000000 && gene.getStop() < 40000000))) {
					hlaGenes.add(gene.getGene());
				}
			}
			LOGGER.info("Excluding " + hlaGenes.size() + " genes");
		} else {
			hlaGenes = null;
		}

		final List<PathwayDatabase> pathwayDatabases = options.getPathwayDatabases();

		final int nrSampleToUseForCorrelation = options.getPermutationGeneCorrelations();
		final int nrSamplesToUseForNullBetas = options.getPermutationPathwayEnrichment();

		final Set<String> nullGwasRuns = genePvaluesNullGwas.getHashCols().keySet();
		if (nullGwasRuns.size() < (nrSampleToUseForCorrelation + nrSamplesToUseForNullBetas)) {
			throw new Exception("Not enough null gwas runs: " + nullGwasRuns.size() + " < " + nrSampleToUseForCorrelation + " + " + nrSamplesToUseForNullBetas);
		}

		Iterator<String> nullGwasRunIterator = nullGwasRuns.iterator();

		final LinkedHashSet<String> sampleToUseForCorrelation = new LinkedHashSet<>(nrSampleToUseForCorrelation);
		for (int i = 0; i < nrSampleToUseForCorrelation; ++i) {
			sampleToUseForCorrelation.add(nullGwasRunIterator.next());
		}

		final LinkedHashSet<String> samplesToUseForNullBetas = new LinkedHashSet<>(nrSamplesToUseForNullBetas);
		for (int i = 0; i < nrSamplesToUseForNullBetas; ++i) {
			samplesToUseForNullBetas.add(nullGwasRunIterator.next());
		}

		final DoubleMatrixDataset<String, String> geneZscoresNullGwasCorrelation = genePvaluesNullGwas.viewColSelection(sampleToUseForCorrelation);
		final DoubleMatrixDataset<String, String> geneZscoresNullGwasNullBetas = genePvaluesNullGwas.viewColSelection(samplesToUseForNullBetas);

		ArrayList<PathwayEnrichments> pathwayEnrichments = new ArrayList<>(pathwayDatabases.size());
		for (PathwayDatabase pathwayDatabase : pathwayDatabases) {
			pathwayEnrichments.add(new PathwayEnrichments(pathwayDatabase, selectedGenes, genes, options.isForceNormalPathwayPvalues(), options.isForceNormalGenePvalues(), genePvalues, geneZscoresNullGwasCorrelation, geneZscoresNullGwasNullBetas, options.getOutputBasePath(), hlaGenes, options.isIgnoreGeneCorrelations(), options.getGenePruningR(), options.getGeneCorrelationWindow(), options.getDebugFolder(), options.getIntermediateFolder()));
		}

		if (options.isSaveOuputAsExcelFiles()) {
			ExcelWriter.saveEnrichmentsToExcel(pathwayEnrichments, options.getOutputBasePath(), genePvalues.getColObjects(), hlaGenes != null);
		} else {
			for (PathwayEnrichments pathwayEnrichment : pathwayEnrichments) {
				//this will make sure z-scores are saved even if make excel is off
				pathwayEnrichment.getEnrichmentZscores();
			}
		}

		LOGGER.info("Completed enrichment analysis for " + pathwayDatabases.size() + " pathway databases");

	}

	private static RandomAccessGenotypeData loadGenotypes(Depict2Options options, List<String> variantsToInclude) throws IOException {
		final RandomAccessGenotypeData referenceGenotypeData;

		final SampleFilter sampleFilter;
		if (options.getGenotypeSamplesFile() != null) {
			sampleFilter = readSampleFile(options.getGenotypeSamplesFile());
		} else {
			sampleFilter = null;
		}

		VariantFilter variantFilter;
		if (variantsToInclude == null) {
			variantFilter = null;
		} else {
			variantFilter = new VariantIdIncludeFilter(new HashSet<>(variantsToInclude));
		}

		if (options.getMafFilter() != 0) {
			VariantFilter mafFilter = new VariantFilterMaf(options.getMafFilter());
			if (variantFilter == null) {
				variantFilter = mafFilter;
			} else {
				variantFilter = new VariantCombinedFilter(variantFilter, mafFilter);
			}
		}

		referenceGenotypeData = options.getGenotypeType().createFilteredGenotypeData(options.getGenotypeBasePath(), 10000, variantFilter, sampleFilter, null, 0.34f);

		return referenceGenotypeData;
	}

	protected static final List<String> readMatrixAnnotations(File file) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(file))).withCSVParser(parser).build();

		ArrayList<String> identifiers = new ArrayList<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			identifiers.add(nextLine[0]);

		}

		return identifiers;

	}

	private static List<Gene> readGenes(File geneFile) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(geneFile))).withCSVParser(parser).withSkipLines(1).build();

		final ArrayList<Gene> genes = new ArrayList<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			genes.add(new Gene(nextLine[0], nextLine[1], Integer.parseInt(nextLine[2]), Integer.parseInt(nextLine[3]), nextLine[5]));

		}

		return genes;

	}

	private static SampleIdIncludeFilter readSampleFile(File sampleFile) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(sampleFile))).withCSVParser(parser).withSkipLines(0).build();

		final HashSet<String> samples = new HashSet<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			samples.add(nextLine[0]);

		}

		return new SampleIdIncludeFilter(samples);

	}

	private static HashSet<String> readVariantFilterFile(File variantFilterFile) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(variantFilterFile))).withCSVParser(parser).withSkipLines(0).build();

		final HashSet<String> variants = new HashSet<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			variants.add(nextLine[0]);

		}

		return variants;

	}

	private static void convertTxtToBin(Depict2Options options) throws IOException, Exception {

		if (options.getConversionColumnIncludeFilter() != null && !options.getConversionColumnIncludeFilter().exists()) {
			throw new FileNotFoundException(options.getConversionColumnIncludeFilter().getAbsolutePath() + " (The system cannot find the file specified)");
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

		if (variantsToExclude.size() > 0) {

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

			matrix = matrix.viewRowSelection(variantsToExclude);

		}

		if (options.getConversionColumnIncludeFilter() != null) {
			List<String> colsToSelect = readMatrixAnnotations(options.getConversionColumnIncludeFilter());
			LOGGER.info("Number of selected columns: " + colsToSelect.size());
			matrix = matrix.viewColSelection(colsToSelect);
		}

		matrix.saveBinary(options.getOutputBasePath());

	}

	private static void convertExpressionMatrixToBin(Depict2Options options) throws IOException, Exception {

		if (options.getConversionColumnIncludeFilter() != null && !options.getConversionColumnIncludeFilter().exists()) {
			throw new FileNotFoundException(options.getConversionColumnIncludeFilter().getAbsolutePath() + " (The system cannot find the file specified)");
		}

		DoubleMatrixDataset<String, String> matrix = DoubleMatrixDataset.loadDoubleTextData(options.getGwasZscoreMatrixPath(), '\t');

		LOGGER.info("Loaded expression matrix with " + matrix.columns() + " samples and " + matrix.rows() + " genes");

		HashSet<String> phenotypesHashSet = new HashSet<>(matrix.columns());

		for (String sample : matrix.getHashCols().keySet()) {
			if (!phenotypesHashSet.add(sample)) {
				throw new Exception("Expression matrix contains a duplicate sample columns: " + sample);
			}
		}

		HashSet<String> variantsHashSet = new HashSet<>(matrix.rows());

		for (String gene : matrix.getHashRows().keySet()) {
			if (!variantsHashSet.add(gene)) {
				throw new Exception("Expression matrix contains a duplicate genes: " + gene);
			}
		}

		if (options.getConversionColumnIncludeFilter() != null) {
			List<String> colsToSelect = readMatrixAnnotations(options.getConversionColumnIncludeFilter());
			LOGGER.info("Number of selected columns: " + colsToSelect.size());
			matrix = matrix.viewColSelection(colsToSelect);
		}

		matrix.normalizeRows();

		LOGGER.info("Normalized genes to have mean 0 and sd 1");

		ArrayList<String> rowNames = matrix.getRowObjects();
		ArrayList<String> nonNanRowNames = new ArrayList<>(matrix.rows());

		rows:
		for (int r = 0; r < matrix.rows(); ++r) {
			for (int c = 0; c < matrix.columns(); ++c) {
				if (Double.isNaN(matrix.getElementQuick(r, c))) {
					continue rows;
				}
			}
			nonNanRowNames.add(rowNames.get(r));

		}

		if (nonNanRowNames.size() < rowNames.size()) {
			matrix = matrix.viewRowSelection(nonNanRowNames);
			LOGGER.info("Removing " + (rowNames.size() - nonNanRowNames.size()) + " rows with NaN after normalizing");
		}

		matrix.saveBinary(options.getOutputBasePath());

	}

	private static void mergeBinMatrix(Depict2Options options) throws IOException, Exception {

		BufferedReader inputReader = new BufferedReader(new FileReader(options.getGwasZscoreMatrixPath()));
		LinkedHashSet<DoubleMatrixDatasetFastSubsetLoader> binMatrices = new LinkedHashSet();
		String line;
		while ((line = inputReader.readLine()) != null) {
			if (line.endsWith(".dat")) {
				line = line.substring(0, line.length() - 4);
			}
			binMatrices.add(new DoubleMatrixDatasetFastSubsetLoader(line));
		}

		LinkedHashSet<String> mergedColNames = new LinkedHashSet(binMatrices.size());
		LinkedHashSet<String> rowNameIntersection = new LinkedHashSet();

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
					while (mergedColNames.contains(newCol + "_" + i++));
					newCol = newCol + "_" + i;
				}
				mergedColNames.add(newCol);
			}

		}

		DoubleMatrixDataset<String, String> mergedData = new DoubleMatrixDataset(rowNameIntersection, mergedColNames);

		int mergedCol = 0;
		for (DoubleMatrixDatasetFastSubsetLoader datasetLoader : binMatrices) {
			DoubleMatrixDataset<String, String> dataset = datasetLoader.loadSubsetOfRowsBinaryDoubleData(rowNameIntersection);

			for (int c = 0; c < dataset.columns(); ++c) {

				mergedData.getCol(mergedCol++).assign(dataset.getCol(c));

			}

		}

		LOGGER.info("Merged data contains: " + mergedData.rows() + " rows and " + mergedData.columns() + " columns");

		mergedData.saveBinary(options.getOutputBasePath());

	}

	private static void mergeConvertTxt(Depict2Options options) throws IOException, Exception {

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

		Map<String, DoubleMatrixDataset<String, String>> summaryStatisticsMap = Collections.synchronizedMap(new HashMap<>());

		gwasSummStats.parallelStream().forEach(summStat -> {

			try {
				final List<String> variantsInZscoreMatrix = DoubleMatrixDataset.readDoubleTextDataRowNames(summStat.getSummStatsFile().getAbsolutePath(), '\t');
				// TODO: Implement support for multiple phenotypes
				//final List<String> phenotypesInZscoreMatrix = DoubleMatrixDataset.readDoubleTextDataColNames(line, '\t');

				//LOGGER.info(variantsInZscoreMatrix.size() + " variants in file: " + fileName);
				HashSet<String> variantsHashSet = new HashSet<>(variantsInZscoreMatrix.size());
				HashSet<String> variantsWithDuplicates = new HashSet<>();

				for (String variant : variantsInZscoreMatrix) {
					if (!variantsHashSet.add(variant)) {
						variantsWithDuplicates.add(variant);
					}
				}

				DoubleMatrixDataset<String, String> matrix;

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
					matrix = summStat.loadSubsetSummStats(variantsHashSet);

				} else {
					matrix = summStat.loadSummStats();
				}

				if (matrix.columns() > 1) {
					throw new Exception("Multi column matrix not supported in merge");
				}

				// Put the matrix in memory
				summaryStatisticsMap.put(summStat.getTrait(), matrix);

				/*			LOGGER.info("Read matrix for " + fileName);
			LOGGER.info("nSNPs: " + matrix.getRowObjects().size());
			LOGGER.info("nCols: " + matrix.getColObjects().size());*/
			} catch (Exception e) {
				throw new RuntimeException(e);
			}
		});

		Set<String> overlappingVariants = new HashSet<>();
		for (DoubleMatrixDataset<String, String> dataset : summaryStatisticsMap.values()) {
			// Put the variant set in memory to avoid having to loop it later on
			if (overlappingVariants.isEmpty()) {
				overlappingVariants.addAll(dataset.getHashRows().keySet());
			} else {
				overlappingVariants.retainAll(dataset.getHashRows().keySet());
			}
		}

		LOGGER.info("Read the following phenotypes: ");

		for (String name : summaryStatisticsMap.keySet()) {
			LOGGER.info(name);
		}

		LOGGER.info("Overlapped variants, retained " + overlappingVariants.size());
		LOGGER.info("Merging over " + summaryStatisticsMap.size() + " phenotypes");

		// Initialize the output matrix
		DoubleMatrixDataset<String, String> finalMergedPvalueMatrix = new DoubleMatrixDataset<>(overlappingVariants, summaryStatisticsMap.keySet());

		int i = 0;
		for (String key : summaryStatisticsMap.keySet()) {
			DoubleMatrixDataset<String, String> curMatrix = summaryStatisticsMap.get(key);
			int j = 0;
			for (String curVariant : overlappingVariants) {
				finalMergedPvalueMatrix.setElementQuick(j, i, curMatrix.viewRow(curVariant).get(0));
				//double curValue= curMatrix.viewRow(curVariant).get(0);
				//finalMergedPvalueMatrix.setElement(curVariant, key, curValue);
				j++;
			}
			i++;
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
						variantsToExclude.add(allVariants.get(c));
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

			for (String dupVariant : variantsToExclude) {
				outputLine[0] = dupVariant;
				excludedVariantWriter.writeNext(outputLine);
			}
			excludedVariantWriter.close();

			LOGGER.info("Encounterd " + variantsToExclude.size() + " variants with NaN values, these are excluded from the conversion. For a full list of excluded variants, see: " + excludedVariantsFile.getPath());

			HashSet<String> variantsToInclude = new HashSet<>(allVariants);
			variantsToInclude.removeAll(variantsToExclude);

			finalMergedPvalueMatrix = finalMergedPvalueMatrix.viewRowSelection(variantsToExclude);

		}

		if (options.getGenotypeBasePath() != null) {
			LOGGER.info("Loading genotype information to convert position  summary statistics to variant IDs");

			File excludedVariantsFile = new File(options.getOutputBasePath() + "_updatedVariantIds.txt");

			final CSVWriter updatedVariantWriter = new CSVWriter(new FileWriter(excludedVariantsFile), '\t', '\0', '\0', "\n");
			final String[] outputLine = new String[2];
			outputLine[0] = "Orginal";
			outputLine[1] = "Updated";
			updatedVariantWriter.writeNext(outputLine);

			RandomAccessGenotypeData genotoypes = loadGenotypes(options, null);

			LinkedHashMap<String, Integer> originalRowHash = finalMergedPvalueMatrix.getHashRows();
			LinkedHashMap<String, Integer> updatedRowHash = new LinkedHashMap(originalRowHash.size());

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

			}

			updatedVariantWriter.close();
			
		}
		
		

		finalMergedPvalueMatrix.saveBinary(options.getOutputBasePath());
	}

	private static void convertBinToTxt(Depict2Options options) throws IOException, Exception {

		DoubleMatrixDataset<String, String> matrix = DoubleMatrixDataset.loadDoubleBinaryData(options.getGwasZscoreMatrixPath());

		String[] columnsToExtract = options.getColumnsToExtract();

		if (columnsToExtract != null) {
			matrix = matrix.viewColSelection(columnsToExtract);
		}

		matrix.save(options.getOutputBasePath() + ".txt");

	}

	private static void convertEqtlToBin(Depict2Options options) throws IOException {

		DoubleMatrixDataset<String, String> matrix = DoubleMatrixDataset.loadTransEqtlExpressionMatrix(options.getGwasZscoreMatrixPath());
		matrix.saveBinary(options.getOutputBasePath());

	}

	public static String formatMsForLog(long ms) {
		//return LOG_TIME_FORMAT.format(new Date(ms));
		return DurationFormatUtils.formatDuration(ms, "H:mm:ss.S");
	}

	private static double[] generateRandomChi2(long numberOfPermutations, int numberOfVariantPerGeneToExpect) {

		final double[] randomChi2;
		if ((numberOfPermutations * numberOfVariantPerGeneToExpect) > Integer.MAX_VALUE - 10) {
			randomChi2 = new double[Integer.MAX_VALUE - 10];
		} else {
			randomChi2 = new double[numberOfVariantPerGeneToExpect * (int) numberOfPermutations];
		}

		final int randomChi2Size = randomChi2.length;
		final int nrThreads = Depict2Options.getNumberOfThreadsToUse();
		final int permPerThread = randomChi2Size / nrThreads;
		final int leftoverPerm = randomChi2Size % nrThreads;

		IntStream.range(0, nrThreads).parallel().forEach(task -> {

			final ThreadLocalRandom rnd = ThreadLocalRandom.current();
			double z;
			for (int p = 0; p < permPerThread; ++p) {
				z = rnd.nextGaussian();
				randomChi2[(task * permPerThread) + p] = z * z;
			}

		});

		if (leftoverPerm > 0) {
			final ThreadLocalRandom rnd = ThreadLocalRandom.current();
			double z;
			for (int p = 0; p < leftoverPerm; ++p) {
				z = rnd.nextGaussian();
				randomChi2[(nrThreads * permPerThread) + p] = z * z;
			}
		}

		return randomChi2;

	}

	private static void correlateGenes(Depict2Options options) throws FileNotFoundException, Exception {

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

			expressionMatrix = DoubleMatrixDataset.loadSubsetOfTextDoubleData(options.getGwasZscoreMatrixPath(), '\t', genes, null);

		} else {
			expressionMatrix = DoubleMatrixDataset.loadDoubleTextData(options.getGwasZscoreMatrixPath(), '\t');
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

	private static void tranposeBinMatrix(Depict2Options options) throws IOException {
		DoubleMatrixDataset<String, String> matrix = DoubleMatrixDataset.loadDoubleBinaryData(options.getGwasZscoreMatrixPath());
		matrix = matrix.viewDice();
		matrix.saveBinary(options.getOutputBasePath());
	}

	private static void doPcaOnBinMatrix(Depict2Options options) throws IOException {

		final DoubleMatrixDataset<String, String> dataset = DoubleMatrixDataset.loadDoubleBinaryData(options.getGwasZscoreMatrixPath());

		PcaColt pcaRes = new PcaColt(dataset, true);

		pcaRes.getEigenvectors().save(options.getOutputBasePath() + "_eigenVectors.txt");
		pcaRes.getEigenValues().save(options.getOutputBasePath() + "_eigenValues.txt");
		pcaRes.getPcs().save(options.getOutputBasePath() + "_pcs.txt");

	}

	protected static class ThreadErrorHandler implements Thread.UncaughtExceptionHandler {

		private final String errorSource;

		public ThreadErrorHandler(String errorSource) {
			this.errorSource = errorSource;
		}

		@Override
		public void uncaughtException(Thread t, Throwable e) {

			System.err.println("Problem running: " + errorSource);
			LOGGER.fatal("Error: " + e.getMessage(), e);
			if (LOGGER.isDebugEnabled()) {
				e.printStackTrace();
			} else {
				System.err.println("See log file for stack trace");
			}
			System.exit(1);
		}
	}

}
