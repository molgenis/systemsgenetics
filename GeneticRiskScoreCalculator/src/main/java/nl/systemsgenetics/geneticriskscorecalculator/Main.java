package nl.systemsgenetics.geneticriskscorecalculator;

import gnu.trove.map.hash.TObjectDoubleHashMap;
import java.io.File;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.regex.Pattern;
import org.apache.commons.cli.ParseException;
import org.apache.log4j.FileAppender;
import org.apache.log4j.Level;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;
import org.apache.log4j.SimpleLayout;
import org.apache.logging.log4j.core.Appender;
import org.apache.logging.log4j.core.LoggerContext;
import org.apache.logging.log4j.core.config.Configurator;
import org.apache.logging.log4j.core.config.LoggerConfig;
import org.apache.logging.log4j.core.layout.PatternLayout;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.multipart.IncompatibleMultiPartGenotypeDataException;

/**
 * Hello world!
 *
 */
public class Main {

	//private static final String VERSION = ResourceBundle.getBundle("verion").getString("application.version");
	private static final String VERSION = "0.0.1";

	private static final String HEADER
			= "  /---------------------------------------\\\n"
			+ "  |     Genetic Risk Score Calculator     |\n"
			+ "  |                                       |\n"
			+ "  |              Rudi Alberts             |\n"
			+ "  |                                       |\n"
			+ "  |             Patrick Deelen            |\n"
			+ "  |                                       |\n"
			+ "  |                                       |\n"
			+ "  | Dasha Zhernakova, Marijke v/d Sijde,  |\n"
			+ "  |   Marc Jan Bonder, Harm-Jan Westra,   |\n"
			+ "  |      Lude Franke, Morris Swertz       |\n"
			+ "  |                                       |\n"
			+ "  |     Genomics Coordination Center      |\n"
			+ "  |        Department of Genetics         |\n"
			+ "  |  University Medical Center Groningen  |\n"
			+ "  \\---------------------------------------/";
	private static final Logger LOGGER = LogManager.getLogger(Main.class);
	private static final DateFormat DATE_TIME_FORMAT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
	private static final Date currentDataTime = new Date();
	private static final Pattern TAB_PATTERN = Pattern.compile("\\t");

	public static void main(String[] args) {

		System.out.println(HEADER);
		System.out.println();
		System.out.println("          --- Version: " + VERSION + " ---");
		System.out.println();
		System.out.println("More information: http://molgenis.org/systemsgenetics");
		System.out.println();

		System.out.println("Current date and time: " + DATE_TIME_FORMAT.format(currentDataTime));
		System.out.println();

		System.out.flush(); //flush to make sure header is before errors
		try {
			Thread.sleep(25); //Allows flush to complete
		} catch (InterruptedException ex) {
		}

		if (args.length == 0) {
			Configuration.printHelp();
			System.exit(1);
		}

		final Configuration configuration;
		try {
			configuration = new Configuration(args);
		} catch (ParseException ex) {
			System.err.println("Invalid command line arguments: ");
			System.err.println(ex.getMessage());
			System.err.println();
			Configuration.printHelp();
			System.exit(1);
			return;
		}

		final File outputFile = configuration.getOutputFile();
		final File outputFolder = outputFile.getParentFile();
		if (outputFolder != null && !outputFolder.exists()) {
			if (!outputFolder.mkdirs()) {
				System.err.println("Failed to create ouput folder at: " + outputFolder.getAbsolutePath());
				System.exit(1);
			}
		}

		final File logFile = new File(outputFile.getAbsolutePath() + ".log");
		startLogging(logFile, true);

		final RandomAccessGenotypeData inputGenotypes;
		try {
			inputGenotypes = configuration.getInputDataType().createGenotypeData(configuration.getInputPaths(), 100);
			System.out.println("Loading reference data complete");
			LOGGER.info("Loading reference data complete");
		} catch (IOException ex) {
			System.err.println("Unable to load reference genotypes file.");
			LOGGER.fatal("Unable to load reference genotypes file.", ex);
			System.exit(1);
			return;
		} catch (IncompatibleMultiPartGenotypeDataException ex) {
			System.err.println("Unable to load reference genotypes file.");
			LOGGER.fatal("Unable to load reference genotypes file.", ex);
			System.exit(1);
			return;
		} catch (GenotypeDataException ex) {
			System.err.println("Unable to load reference genotypes file.");
			LOGGER.fatal("Unable to load reference genotypes file.", ex);
			System.exit(1);
			return;
		}
		int debug = 0;

		int index = 0;
		int index2 = 0;
		for (String mySample : inputGenotypes.getSampleNames()) {
			if (index < 5) {
				System.out.println("Samplename: " + mySample);
			}
			index++;
		}

		final String onlyCount = configuration.getOnlyCount();
		final String harmonizedData = configuration.getHarmonizedData();
		final double inclusionThreshold = configuration.getInclusionThreshold();

		GwasCatalogLoader gwasCatalogLoader = new GwasCatalogLoader();
		List<GeneticRiskScoreCalculator> geneticRiskScoreCalculators = gwasCatalogLoader.getGeneticRiskScoreCalculators(configuration.getRisksnpsFile().toString());

		List<String> phenotypes = new ArrayList<String>();
		for (GeneticRiskScoreCalculator calculator : geneticRiskScoreCalculators) {
			phenotypes.add(calculator.getPhenotype());
			System.out.println("USING PHENOTYPE: **" + calculator.getPhenotype() + "**");
		}

		// List<String> phenotypes = Arrays.asList("Height", "Migraine"); // fixed size list
		RiskScoreMatrix riskScoreMatrix = new RiskScoreMatrix(inputGenotypes.getSampleNames(), phenotypes);

		for (GeneticRiskScoreCalculator calculator : geneticRiskScoreCalculators) {
			TObjectDoubleHashMap<String> riskScores = calculator.calculateRiskScores(inputGenotypes, onlyCount, harmonizedData, inclusionThreshold);
			index = 0;
			for (String sample : inputGenotypes.getSampleNames()) {
				riskScoreMatrix.setRiskScore(sample, calculator.getPhenotype(), riskScores.get(sample));
				if (index < 3) {
					System.out.println("sample:" + sample);
					System.out.println("phenotype:" + calculator.getPhenotype());
					System.out.println("score:" + riskScores.get(sample));
				}
				index++;
				//System.out.println("sample:" + calculator.getPhenotype());
			}

		}

		//  System.out.println("GET A SCORE: " + riskScoreMatrix.getRiskScore("LLDeep_1094", "Height"));
		//  System.out.println("GET A SCORE: " + riskScoreMatrix.getRiskScore("LLDeep_0727", "Height"));
		//  System.out.println("GET A SCORE: " + riskScoreMatrix.getRiskScore("LLDeep_1094", "Migraine"));
		//  System.out.println("GET A SCORE: " + riskScoreMatrix.getRiskScore("LLDeep_0727", "Migraine"));
		System.out.println("GET A SCORE ROWS: " + riskScoreMatrix.rows());
		System.out.println("GET A SCORE COLS: " + riskScoreMatrix.cols());

		try {
			riskScoreMatrix.save(outputFile);
		} catch (IOException ex) {
			System.err.println("Could not save output file to: " + outputFile.getAbsolutePath());
			LOGGER.fatal("Could not save output file to: " + outputFile.getAbsolutePath(), ex);
			System.exit(1);
			return;
		}
		System.out.println("Risk score calculation complete");
		LOGGER.info("Risk score calculation complete");

	}
// get sample dosages
	// 0,1,2   ....    1.5   1.4   1.33434   
	// imputation not sure .... in between aa ab bb

	private static void startLogging(File logFile, boolean debugMode) {

		try {

			if (debugMode) {
				Configurator.setRootLevel(org.apache.logging.log4j.Level.DEBUG);
			} else {
				Configurator.setRootLevel(org.apache.logging.log4j.Level.INFO);
			}
			LoggerContext context = LoggerContext.getContext(false);
			org.apache.logging.log4j.core.config.Configuration config = context.getConfiguration();

			PatternLayout loggingLayoutFull = PatternLayout.newBuilder()
					.withPattern("[%level] %d{ABSOLUTE} - %c{1} - %msg%n")
					.build();

			PatternLayout loggingLayoutReduced = PatternLayout.newBuilder()
					.withPattern("%msg%n")
					.build();

			// Log file appender
			org.apache.logging.log4j.core.appender.FileAppender file = org.apache.logging.log4j.core.appender.FileAppender.newBuilder()
					.setName("file")
					.setLayout(loggingLayoutFull)
					.withFileName(logFile.getCanonicalPath())
					.build();
			file.start();

			// Make sure any existing loggers are removed
			for (Appender appender : context.getRootLogger().getAppenders().values()) {
				context.getRootLogger().removeAppender(appender);
			}

			// Add the appenders to the root logger
			Logger rootLogger = context.getRootLogger();
			LoggerConfig rootLoggerConfig = config.getRootLogger();

			rootLoggerConfig.addAppender(file, org.apache.logging.log4j.Level.INFO, null);
			config.addLogger(rootLogger.getName(), rootLoggerConfig);

			context.updateLoggers(config);

		} catch (IOException e) {
			System.err.println("Failed to create logger: " + e.getMessage());
			System.exit(1);
		}

		LOGGER.info(
				"\n" + HEADER);
		LOGGER.info(
				"Version: " + VERSION);
		LOGGER.info(
				"Current date and time: " + DATE_TIME_FORMAT.format(currentDataTime));
		LOGGER.info(
				"Log level: " + LOGGER.getLevel());

		System.out.println(
				"Started logging");
		System.out.println();
	}

}
