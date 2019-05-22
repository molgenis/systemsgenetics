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
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.ResourceBundle;
import java.util.concurrent.ThreadLocalRandom;
import java.util.stream.IntStream;
import nl.systemsgenetics.depict2.development.ExtractCol;
import nl.systemsgenetics.depict2.development.First1000qtl;
import org.apache.commons.cli.ParseException;
import org.apache.commons.lang.time.DurationFormatUtils;
import org.apache.log4j.ConsoleAppender;
import org.apache.log4j.FileAppender;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.SimpleLayout;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.multipart.IncompatibleMultiPartGenotypeDataException;
import org.molgenis.genotype.sampleFilter.SampleFilter;
import org.molgenis.genotype.sampleFilter.SampleIdIncludeFilter;
import org.molgenis.genotype.tabix.TabixFileNotFoundException;
import org.molgenis.genotype.variantFilter.VariantIdIncludeFilter;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.ZScores;

/**
 *
 * @author Patrick Deelen
 */
public class Depict2 {

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
			FileAppender logFileAppender = new FileAppender(new SimpleLayout(), options.getLogFile().getCanonicalPath(), options.getMode() == Depict2Mode.RUN2 || options.getMode() == Depict2Mode.RUN3);
			ConsoleAppender logConsoleInfoAppender = new ConsoleAppender(new InfoOnlyLogLayout());
			Logger.getRootLogger().removeAllAppenders();
			Logger.getRootLogger().addAppender(logFileAppender);

			LOGGER.info("DEPICT" + VERSION);
			LOGGER.info("Current date and time: " + startDateTime);

			Logger.getRootLogger().addAppender(logConsoleInfoAppender);

			if (options.isDebugMode()) {
				Logger.getRootLogger().setLevel(Level.DEBUG);
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
				case FIRST1000:
					First1000qtl.printFirst1000(options);
					break;
				case RUN:
					run(options);
					break;
				case RUN2:
					run2(options, null, null, null, null);
					break;
				case RUN3:
					run3(options, null, null, null, null);
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

		LOGGER.info("Number of phenotypes in GWAS matrix: " + phenotypesInZscoreMatrix.size());
		LOGGER.info("Number of variants in GWAS matrix: " + variantsInZscoreMatrix.size());

		RandomAccessGenotypeData referenceGenotypeData = loadGenotypes(options, variantsInZscoreMatrix);

		LOGGER.info("Done loading genotype data");

		List<Gene> genes = readGenes(options.getGeneInfoFile());

		LOGGER.info("Loaded " + genes.size() + " genes");

		double[] randomChi2 = generateRandomChi2(options.getNumberOfPermutations(), 500);

		LOGGER.info("Prepared reference null distribution with " + randomChi2.length + " values");

		GenePvalueCalculator gpc = new GenePvalueCalculator(options.getGwasZscoreMatrixPath(), referenceGenotypeData, genes, options.getWindowExtend(), options.getMaxRBetweenVariants(), options.getNumberOfPermutations(), options.getOutputBasePath(), randomChi2, options.correctForLambdaInflation());

		DoubleMatrixDataset<String, String> genePvalues = gpc.getGenePvalues();
		DoubleMatrixDataset<String, String> genePvaluesNullGwas = gpc.getGenePvaluesNullGwas();
		DoubleMatrixDataset<String, String> geneVariantCount = gpc.getGeneVariantCount();

		genePvalues.save(options.getOutputBasePath() + "_genePvalues.txt");
		genePvaluesNullGwas.save(options.getOutputBasePath() + "_genePvaluesNullGwas.txt");
		geneVariantCount.save(options.getOutputBasePath() + "_geneVariantCount.txt");
		gpc.getGeneMaxPermutationCount().save(options.getOutputBasePath() + "_geneMaxPermutationUsed.txt");
		gpc.getGeneRuntime().save(options.getOutputBasePath() + "_geneRuntime.txt");

		LOGGER.info("Gene p-values saved. If needed the analysis can be resummed from this point using --mode RUN2 and exactly the same output path and genes file");

		run2(options, genePvalues, genePvaluesNullGwas, geneVariantCount, genes);

	}

	/**
	 * Part2 of depict. Matrices may be null and then they will be loaded form
	 * output path. If matrices is null genes will be reloaded
	 *
	 * All matrices will have genes on the rows in the same order
	 *
	 * @param options
	 * @param genePvalues
	 * @param genePvaluesNullGwas
	 * @param geneVariantCount
	 * @throws IOException
	 * @throws Exception
	 */
	private static void run2(Depict2Options options, DoubleMatrixDataset<String, String> genePvalues, DoubleMatrixDataset<String, String> genePvaluesNullGwas, DoubleMatrixDataset<String, String> geneVariantCount, List<Gene> genes) throws IOException, Exception {

		if (options.getMode() == Depict2Mode.RUN2) {
			LOGGER.info("Continuing previous analysis by loading gene p-values");
			genePvalues = DoubleMatrixDataset.loadDoubleTextData(options.getOutputBasePath() + "_genePvalues.txt", '\t');
			genePvaluesNullGwas = DoubleMatrixDataset.loadDoubleTextData(options.getOutputBasePath() + "_genePvaluesNullGwas.txt", '\t');
			geneVariantCount = DoubleMatrixDataset.loadDoubleTextData(options.getOutputBasePath() + "_geneVariantCount.txt", '\t');
			LOGGER.info("Gene p-values loaded");
			genes = readGenes(options.getGeneInfoFile());
			LOGGER.info("Loaded " + genes.size() + " genes");
		}

		//Identify genes with atleast one variant in window
		final LinkedHashSet<String> selectedGenes = new LinkedHashSet<>();
		final ArrayList<String> allGenes = geneVariantCount.getRowObjects();
		final int totalGeneCount = allGenes.size();
		for (int g = 0; g < totalGeneCount; ++g) {
			if (geneVariantCount.getElementQuick(g, 0) > 0) {
				selectedGenes.add(allGenes.get(g));
			}
		}
		LOGGER.info("Number of genes with atleast one variant in specified window: " + selectedGenes.size());

		//Exclude genes with no variants
		genePvalues = genePvalues.viewRowSelection(selectedGenes);
		genePvaluesNullGwas = genePvaluesNullGwas.viewRowSelection(selectedGenes);
		//geneVariantCount = geneVariantCount.viewRowSelection(selectedGenes);

		//Gene weight will have same order as other matrices
		List<DoubleMatrixDataset<String, String>> invCorMatrixPerChrArm = CalculateGeneInvCorMatrix.CalculateGeneInvCorMatrix(genePvaluesNullGwas, genes, options);

		

		if (options.getPathwayDatabases().isEmpty()) {
			LOGGER.info("Gene weights saved. The analysis will now stop since no pathway databases are provided. Use --mode RUN3 and exactly the same output path and genes file to continue");
		} else {
			LOGGER.info("Gene weights saved. If needed the analysis can be resummed from this point using --mode RUN3 and exactly the same output path and genes file");
			run3(options, genePvalues, genePvaluesNullGwas, genes, invCorMatrixPerChrArm);
		}

	}

	private static void run3(Depict2Options options, DoubleMatrixDataset<String, String> genePvalues, DoubleMatrixDataset<String, String> genePvaluesNullGwas, List<Gene> genes, List<DoubleMatrixDataset<String, String>> invCorMatrixPerChrArm) throws IOException, Exception {

		if (options.getMode() == Depict2Mode.RUN3) {
			LOGGER.info("Continuing previous analysis by loading gene p-values and gene weigthts");
			genePvalues = DoubleMatrixDataset.loadDoubleTextData(options.getOutputBasePath() + "_genePvalues.txt", '\t');
			genePvaluesNullGwas = DoubleMatrixDataset.loadDoubleTextData(options.getOutputBasePath() + "_genePvaluesNullGwas.txt", '\t');
			LOGGER.info("Gene p-values loaded");
			genes = readGenes(options.getGeneInfoFile());
			LOGGER.info("Loaded " + genes.size() + " genes");
			
			final Map<String, ArrayList<String>> chrArmToGeneMapping = createChrArmGeneMapping(genes, genePvaluesNullGwas.getHashRows());
			
			geneInvCorMatrix = DoubleMatrixDataset.loadDoubleTextData(options.getOutputBasePath() + "_GeneCorMatrix.txt", '\t');
			LOGGER.info("Gene weights loaded");

			//In run2 we only calculate weigths for genes with atleast one variant in the GWAS. We have to redo this selection
			genePvalues = genePvalues.viewRowSelection(geneInvCorMatrix.getHashRows().keySet());
			genePvaluesNullGwas = genePvaluesNullGwas.viewRowSelection(geneInvCorMatrix.getHashRows().keySet());

		}

		List<PathwayDatabase> pathwayDatabases = options.getPathwayDatabases();

		DoubleMatrix2D matrix = genePvalues.getMatrix();

		//Inplace convert gene p-values to z-scores
		for (int r = 0; r < matrix.rows(); ++r) {
			for (int c = 0; c < matrix.columns(); ++c) {
				matrix.setQuick(r, c, -ZScores.pToZTwoTailed(matrix.getQuick(r, c)));
			}
		}

//		genePvalues = genePvalues.viewDice().createRowForceNormalDuplicate().viewDice();
		DoubleMatrix2D matrixNull = genePvaluesNullGwas.getMatrix();

		for (int r = 0; r < matrixNull.rows(); ++r) {
			for (int c = 0; c < matrixNull.columns(); ++c) {
				matrixNull.setQuick(r, c, -ZScores.pToZTwoTailed(matrixNull.getQuick(r, c)));
			}
		}

//		genePvaluesNullGwas = genePvaluesNullGwas.viewDice().createRowForceNormalDuplicate().viewDice();
		HashMap<PathwayDatabase, DoubleMatrixDataset<String, String>> enrichments = PathwayEnrichments.performEnrichmentAnalysis(genePvalues, genePvaluesNullGwas, geneInvCorMatrix, pathwayDatabases, options.getOutputBasePath(), null);

		PathwayEnrichments.saveEnrichmentsToExcel(pathwayDatabases, options.getOutputBasePath(), enrichments, genePvalues.getColObjects(), false);

		LOGGER.info("Completed enrichment analysis for " + pathwayDatabases.size() + " pathway databases");

		HashSet<String> hlaGenes = new HashSet<>();
		for (Gene gene : genes) {
			if (gene.getChr().equals("6") && ((gene.getStart() > 20000000 && gene.getStart() < 40000000) || (gene.getStop() > 20000000 && gene.getStop() < 40000000))) {
				hlaGenes.add(gene.getGene());
			}
		}

		enrichments = PathwayEnrichments.performEnrichmentAnalysis(genePvalues, genePvaluesNullGwas, geneInvCorMatrix, pathwayDatabases, options.getOutputBasePath(), hlaGenes);

		PathwayEnrichments.saveEnrichmentsToExcel(pathwayDatabases, options.getOutputBasePath(), enrichments, genePvalues.getColObjects(), true);
		
		LOGGER.info("Completed enrichment without " + hlaGenes.size() + " gene in HLA region for " + pathwayDatabases.size() + " pathway databases");

	}

	private static RandomAccessGenotypeData loadGenotypes(Depict2Options options, List<String> variantsToInclude) throws IOException {
		final RandomAccessGenotypeData referenceGenotypeData;

		final SampleFilter sampleFilter;
		if (options.getGenotypeSamplesFile() != null) {
			sampleFilter = readSampleFile(options.getGenotypeSamplesFile());
		} else {
			sampleFilter = null;
		}

		referenceGenotypeData = options.getGenotypeType().createFilteredGenotypeData(options.getGenotypeBasePath(), 10000, new VariantIdIncludeFilter(new HashSet<>(variantsToInclude)), sampleFilter, null, 0.34f);

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

			File excludedVariantsFile = new File(options.getOutputBasePath() + "_excludedVariants.txt");

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

		if (options.isPvalueToZscore()) {
			DoubleMatrix2D matrixContent = matrix.getMatrix();

			int rows = matrixContent.rows();
			int cols = matrixContent.columns();

			for (int r = 0; r < rows; ++r) {
				for (int c = 0; c < cols; ++c) {

					matrixContent.setQuick(r, c, ZScores.pToZTwoTailed(matrixContent.getQuick(r, c)));

				}
			}

		}

		if (options.getConversionColumnIncludeFilter() != null) {
			List<String> colsToSelect = readMatrixAnnotations(options.getConversionColumnIncludeFilter());
			LOGGER.info("Number of selected columns: " + colsToSelect.size());
			matrix = matrix.viewColSelection(colsToSelect);
		}

		matrix.saveBinary(options.getOutputBasePath());

	}
	
	private static void convertBinToTxt(Depict2Options options) throws IOException, Exception {
		
		DoubleMatrixDataset<String, String> matrix = DoubleMatrixDataset.loadDoubleBinaryData(options.getGwasZscoreMatrixPath());
		matrix.save(options.getOutputBasePath());
		
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
