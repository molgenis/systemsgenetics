package nl.systemsgenetics.downstreamer;

import nl.systemsgenetics.downstreamer.runners.PreparePermutations;
import nl.systemsgenetics.downstreamer.io.ExcelWriter;
import nl.systemsgenetics.downstreamer.pathway.PredictedPathwayAnnotations;
import nl.systemsgenetics.downstreamer.runners.options.*;
import nl.systemsgenetics.downstreamer.runners.*;
import org.apache.commons.cli.ParseException;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.lang.time.DurationFormatUtils;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.Appender;
import org.apache.logging.log4j.core.LoggerContext;
import org.apache.logging.log4j.core.appender.ConsoleAppender;
import org.apache.logging.log4j.core.appender.FileAppender;
import org.apache.logging.log4j.core.config.Configuration;
import org.apache.logging.log4j.core.config.Configurator;
import org.apache.logging.log4j.core.config.LoggerConfig;
import org.apache.logging.log4j.core.layout.PatternLayout;


import java.io.File;
import java.io.IOException;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.ResourceBundle;

/**
 * Entry point and error handler for Downstreamer
 *
 * @author Patrick Deelen
 */
public class Downstreamer {

	public static final DecimalFormat LARGE_INT_FORMAT = new DecimalFormat("###,###");
	public static final String VERSION = ResourceBundle.getBundle("verion").getString("application.version");
	public static final DateFormat DATE_TIME_FORMAT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
	public static final DateFormat TIME_FORMAT = new SimpleDateFormat("HH:mm:ss");

	private static final Logger LOGGER = LogManager.getLogger(Downstreamer.class);
	private static final String HEADER
			= "  /---------------------------------------\\\n"
			+ "  |             Downstreamer              |\n"
			+ "  |                                       |\n"
			+ "  |  University Medical Center Groningen  |\n"
			+ "  \\---------------------------------------/";

	/**
	 * Main function for Downstreamer analysis, specific tasks are implemented
	 * in respective runner classes.
	 *
	 * @param args the command line arguments
	 * @throws InterruptedException
	 */
	public static void main(String[] args) throws InterruptedException {
		System.out.println();
		System.out.println(HEADER);
		System.out.println();
		System.out.println("       --- Version: " + VERSION + " ---");
		System.out.println();
		System.out.println("More information: https://github.com/molgenis/systemsgenetics/wiki/Downstreamer");
		System.out.println();

		Date currentDataTime = new Date();
		String startDateTime = DATE_TIME_FORMAT.format(currentDataTime);
		System.out.println("Current date and time: " + startDateTime);
		System.out.println();

		System.out.flush(); //flush to make sure header is before errors
		Thread.sleep(25); //Allows flush to complete

		final OptionsBase options;

		// If no arguments are provided print help
		if (args.length == 0 || args[0].equals("-h") && args.length == 1 || args[0].equals("--help") && args.length == 1) {
			OptionsBase.printHelp();
			return;
		}

		// Parse general options
		try {
			options = new OptionsBase(args);
		} catch (ParseException ex) {
			System.err.println("Error parsing commandline: " + ex.getMessage());
			OptionsBase.printHelp();
			return;
		}

		
		
		// If the output folder does not exist, create it
		if (options.getLogFile().getParentFile() != null && !options.getLogFile().getParentFile().isDirectory()) {
			if (!options.getLogFile().getParentFile().mkdirs()) {
				System.err.println("Failed to create output folder: " + options.getLogFile().getParent());
				System.exit(1);
			}
		}

		// Initialize logfile
		try {

			Level loggingLevel = Level.INFO;
			if (options.isDebugMode()) {
				loggingLevel = Level.DEBUG;
				OptionsBase.getDebugFolder().mkdir();
			}

			initializeLoggers(loggingLevel, options.getLogFile(), startDateTime);

		} catch (IOException e) {
			System.err.println("Failed to create logger: " + e.getMessage());
			System.exit(1);
		}

		// Run mode
		// This allows for dynamic help messages so there is not a full dump of all possible options
		// which I found really annoying when trying to run DS myself.
		HelpProvider helpProvider = OptionsBase::printHelp;

		// For legacy modes
		DownstreamerOptionsDeprecated optionsDeprecated;
		try {

			// Set the function that provides the help text
			switch (options.getMode()) {
				case REGRESS:
					helpProvider = OptionsModeRegress::printHelp;
					break;
				case ENRICH:
					helpProvider = OptionsModeEnrichment::printHelp;
					break;
				case PREPARE_PERMUTATIONS:
					helpProvider = OptionsModePerparePermutations::printHelp;
					break;
				case CONVERT_BIN:
				case CONVERT_TXT:
				case CONVERT_EQTL:
				case CONVERT_GTEX:
				case CONVERT_PTOZSCORE:
				case CONVERT_MERGE_BIN:
				case CONVERT_DAT_TO_DATG:
				case CONVERT_EXP:
				case PREPARE_GWAS:
				case PREPARE_GWAS_MERGE:
				case MATRIX_SUBSET:
				case MATRIX_TRANSPOSE:
					helpProvider = OptionsModeConvert::printHelp;
					break;
				case COREG_CORRELATE_GENES:
				case COREG_RTOZSCORE:
				case COREG_PCA:
				case COREG_REMOVE_CIS_COEXP:
				case COREG_INVESTIGATE_NETWORK:
					helpProvider = OptionsModeCoreg::printHelp;
					break;
				// Legacy options to update
				case CREATE_EXCEL:
				case PRIO_GENE_ENRICH:
				case EXPAND_PATHWAYS:
				case STEP1:
				case STEP2:
					helpProvider = DownstreamerOptionsDeprecated::printHelp;
					break;
				default:
					helpProvider = OptionsBase::printHelp;
					break;
			}

			// Execute the right mode
			switch (options.getMode()) {
				// New main analyses
				case REGRESS:
					DownstreamerRegressionEngine.run(new OptionsModeRegress(args));
					break;
				case ENRICH:
					OptionsModeEnrichment optionsRegress = new OptionsModeEnrichment(args);
					DownstreamerStep2Results results = DownstreamerEnrichment.enrichmentAnalysis(optionsRegress);
					ExcelWriter writer = new ExcelWriter(results.getGenePvalues().getColObjects(), optionsRegress);
					writer.saveStep2Excel(results);
					break;

				case PREPARE_PERMUTATIONS:
					PreparePermutations.run(new OptionsModePerparePermutations(args));
					break;

				// Converter utils that share OptionsModeConvert()
				case CONVERT_TXT:
					options.printOptions();
					DownstreamerConverters.convertTxtToBin(new OptionsModeConvert(args));
					break;
				case CONVERT_BIN:
					options.printOptions();
					DownstreamerConverters.convertBinToTxt(new OptionsModeConvert(args));
					break;
				case CONVERT_MERGE_BIN:
					options.printOptions();
					DownstreamerConverters.mergeBinMatrix(new OptionsModeConvert(args));
					break;
				case CONVERT_EXP:
					options.printOptions();
					DownstreamerConverters.convertExpressionMatrixToBin(new OptionsModeConvert(args));
					break;
				case CONVERT_DAT_TO_DATG:
					options.printOptions();
					DownstreamerConverters.convertDatToDatg(new OptionsModeConvert(args));
					break;
				case CONVERT_GTEX:
					options.printOptions();
					DownstreamerConverters.convertGct(new OptionsModeConvert(args));
					break;
				case CONVERT_EQTL:
					options.printOptions();
					DownstreamerConverters.convertEqtlToBin(new OptionsModeConvert(args));
					break;
				case CONVERT_PTOZSCORE:
					options.printOptions();
					DownstreamerConverters.convertPvalueToZscore(new OptionsModeConvert(args));
					break;
				case PREPARE_GWAS:
					options.printOptions();
					DownstreamerConverters.prepareGwasSummaryStatistics(new OptionsModeConvert(args));
					break;
				case PREPARE_GWAS_MERGE:
					options.printOptions();
					DownstreamerConverters.prepareGwasSummaryStatisticsMerge(new OptionsModeConvert(args));
					break;
				case MATRIX_SUBSET:
					options.printOptions();
					DownstreamerConverters.subsetMatrix(new OptionsModeConvert(args));
					break;
				case MATRIX_TRANSPOSE:
					options.printOptions();
					DownstreamerConverters.tranposeBinMatrix(new OptionsModeConvert(args));
					break;

				// Utilities related to making the co-regulation network. Share OptionsModeCoreg()
				case COREG_CORRELATE_GENES:
					options.printOptions();
					CoregulationUtilities.coregCorrelateGenes(new OptionsModeCoreg(args));
					break;
				case COREG_RTOZSCORE:
					options.printOptions();
					CoregulationUtilities.coregConvertRtoZscore(new OptionsModeCoreg(args));
					break;
				case COREG_PCA:
					options.printOptions();
					CoregulationUtilities.coregDoPcaOnBinMatrix(new OptionsModeCoreg(args));
					break;
				case COREG_REMOVE_CIS_COEXP:
					options.printOptions();
					CoregulationUtilities.coregRemoveLocalGeneCorrelations(new OptionsModeCoreg(args));
					break;
				case COREG_INVESTIGATE_NETWORK:
					options.printOptions();
					CoregulationUtilities.coregInvestigateNetwork(new OptionsModeCoreg(args));
					break;

				// Testing utils
				case TEST_DECOMP:
					options.printOptions();
					DownstreamerUtilities.testEigenDecomposition(new OptionsTesting(args));
					break;

				// Legacy modes, These will be either replaced or substantially updated
				case CREATE_EXCEL:
					options.printOptions();
					DownstreamerUtilities.generateExcelFromIntermediates(new DownstreamerOptionsDeprecated(args));
					break;
				case PRIO_GENE_ENRICH:
					options.printOptions();
					PathwayDatabaseEnrichments.testPredictionPerformance(new DownstreamerOptionsDeprecated(args));
					break;
				case EXPAND_PATHWAYS:
					options.printOptions();
					PredictedPathwayAnnotations.expandAnnotations(new DownstreamerOptionsDeprecated(args));
					break;
				case STEP1:
//					optionsDeprecated = new DownstreamerOptionsDeprecated(args);
//					final DownstreamerStep1Results step1Res = DownstreamerMainAnalysis.step1(optionsDeprecated);
//					final DownstreamerStep2Results step2Res;
//
//					if (optionsDeprecated.getPathwayDatabases().isEmpty()) {
//						LOGGER.info("The analysis will now stop since no pathway databases are provided. Use --mode STEP2 and exactly the same output path and genes file to continue");
//						break;
//					} else {
//						step2Res = DownstreamerMainAnalysis.step2(optionsDeprecated, step1Res);
//					}
//
//					if (step2Res != null) {
//						ExcelWriter writer = null;
//						if (optionsDeprecated.isSaveOuputAsExcelFiles()) {
//							writer = new ExcelWriter(step2Res.getGenePvalues().getColObjects(), optionsDeprecated);
//							writer.saveStep2Excel(step2Res);
//						}
//
//						if (optionsDeprecated.isAssignPathwayGenesToCisWindow()) {
//							final DownstreamerStep3Results step3Res = DownstreamerMainAnalysis.step3(optionsDeprecated);
//							if (writer != null) {
//								writer.saveStep3Excel(step2Res, step3Res);
//							}
//						}
//					}
//
//					break;

				case STEP2:
//					optionsDeprecated = new DownstreamerOptionsDeprecated(args);
//					step2Res = DownstreamerMainAnalysis.step2(optionsDeprecated, null);
//
//					ExcelWriter writer = null;
//					if (optionsDeprecated.isSaveOuputAsExcelFiles()) {
//						writer = new ExcelWriter(step2Res.getGenePvalues().getColObjects(), optionsDeprecated);
//						writer.saveStep2Excel(step2Res);
//					}
//					if (optionsDeprecated.isAssignPathwayGenesToCisWindow()) {
//						final DownstreamerStep3Results step3Res = DownstreamerMainAnalysis.step3(optionsDeprecated);
//						if (writer != null) {
//							writer.saveStep3Excel(step2Res, step3Res);
//						}
//					}
//					break;
				default:
					throw new Exception("Mode not yet converted or is deprecated.");
			}
		} catch (ParseException e) {
			System.err.println("Problem running mode: " + options.getMode());
			LOGGER.fatal("Error: ", e);
			helpProvider.printHelp();
			System.exit(1);
		} catch (Exception e) {
			System.err.println("Problem running mode: " + options.getMode());

			LOGGER.fatal("Error: ", e);
			if (LOGGER.isDebugEnabled()) {
				e.printStackTrace();
			} else {
				System.err.println("See log file for stack trace");
			}

			System.exit(1);
		}
		LOGGER.info("Analysis completed");

		currentDataTime = new Date();
		LOGGER.info("Current date and time: " + DATE_TIME_FORMAT.format(currentDataTime));
		LOGGER.info(StringUtils.repeat("-", 80));
		LOGGER.info("");
	}

	public static void initializeLoggers(Level loggingLevel, File logFile, String startDateTime) throws IOException {
		Configurator.setRootLevel(loggingLevel);
		org.apache.logging.log4j.core.LoggerContext context = LoggerContext.getContext(false);
		Configuration config = context.getConfiguration();

		PatternLayout loggingLayoutFull = PatternLayout.newBuilder()
				.withPattern("[%level]\t%d{ABSOLUTE} - %msg%n")
				.build();

		PatternLayout loggingLayoutReduced = PatternLayout.newBuilder()
				.withPattern("%msg %throwable{short.message}%n")
				.build();

		// Stdout appender
		ConsoleAppender stdOut = ConsoleAppender.newBuilder()
				.setName("stdout")
				.setLayout(loggingLayoutReduced)
				.build();
		stdOut.start();

		// Log file appender
		FileAppender file = FileAppender.newBuilder()
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

		rootLoggerConfig.addAppender(file, loggingLevel, null);
		config.addLogger(rootLogger.getName(), rootLoggerConfig);

		context.updateLoggers(config);

		// Create a seperator in the log, so seperate runs are easier to spot
		LOGGER.info("");
		LOGGER.info("Downstreamer " + VERSION);
		LOGGER.info("Current date and time: " + startDateTime);

		rootLoggerConfig.addAppender(stdOut, Level.INFO, null);
	}

	public static double memUsage() {
		Runtime rt = Runtime.getRuntime();
		double used = (rt.totalMemory() - rt.freeMemory()) / 1024.0 / 1024 / 1024;
		return used;
	}

	public static void logDebugMem(String message, Logger logger) {
		logger.debug("[mem: " + String.format("%,.2f", memUsage()) + "G] " + message);
	}

	public static String formatMsForLog(long ms) {
		//return LOG_TIME_FORMAT.format(new Date(ms));
		return DurationFormatUtils.formatDuration(ms, "H:mm:ss.S");
	}

	public static class ThreadErrorHandler implements Thread.UncaughtExceptionHandler {

		private final String errorSource;

		public ThreadErrorHandler(String errorSource) {
			this.errorSource = errorSource;
		}

		@Override
		public void uncaughtException(Thread t, Throwable e) {

			System.err.println("Problem running: " + errorSource);
			LOGGER.fatal("Error: ", e);
			if (LOGGER.isDebugEnabled()) {
				e.printStackTrace();
			} else {
				System.err.println("See log file for stack trace");
			}
			System.exit(1);
		}
	}

}
