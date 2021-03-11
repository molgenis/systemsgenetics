package nl.systemsgenetics.downstreamer;

import nl.systemsgenetics.downstreamer.runners.TestCoregulationPerformance;
import nl.systemsgenetics.downstreamer.runners.NetworkProperties;
import nl.systemsgenetics.downstreamer.runners.DownstreamerUtilities;
import nl.systemsgenetics.downstreamer.runners.DownstreamerMainAnalysis;
import nl.systemsgenetics.downstreamer.runners.PruneToIndependentTopHits;
import nl.systemsgenetics.downstreamer.runners.DownstreamerConverters;
import java.io.File;
import java.io.IOException;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.*;

import nl.systemsgenetics.downstreamer.development.CorrelateExpressionToPredictions;
import nl.systemsgenetics.downstreamer.development.First1000qtl;
import nl.systemsgenetics.downstreamer.io.ExcelWriter;

import org.apache.commons.cli.ParseException;
import org.apache.commons.lang.time.DurationFormatUtils;
import org.apache.log4j.ConsoleAppender;
import org.apache.log4j.FileAppender;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.SimpleLayout;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.multipart.IncompatibleMultiPartGenotypeDataException;
import org.molgenis.genotype.tabix.TabixFileNotFoundException;

/**
 * Entry point and error handler for Downstreamer
 *
 * @author Patrick Deelen
 */
public class Downstreamer {

	public static final DecimalFormat LARGE_INT_FORMAT = new DecimalFormat("###,###");
	public static final String VERSION = ResourceBundle.getBundle("verion").getString("application.version");
	private static final DateFormat DATE_TIME_FORMAT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
	private static final Logger LOGGER = Logger.getLogger(Downstreamer.class);
	private static final String HEADER
			= "  /---------------------------------------\\\n"
			+ "  |             Downstreamer              |\n"
			+ "  |                                       |\n"
			+ "  |  University Medical Center Groningen  |\n"
			+ "  \\---------------------------------------/";

	/**
	 * Main function for Downstreamer analysis, specific tasks are implemented in
	 * respective runner classes.
	 *
	 * @param args the command line arguments
	 * @throws java.lang.InterruptedException
	 */
	public static void main(String[] args) throws InterruptedException {

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

		DownstreamerOptions options;

		if (args.length == 0) {
			DownstreamerOptions.printHelp();
			return;
		}

		try {
			options = new DownstreamerOptions(args);
		} catch (ParseException ex) {
			System.err.println("Error parsing commandline: " + ex.getMessage());
			DownstreamerOptions.printHelp();
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
			FileAppender logFileAppender = new FileAppender(new SimpleLayout(), options.getLogFile().getCanonicalPath(), options.getMode() != DownstreamerMode.STEP1);
			ConsoleAppender logConsoleInfoAppender = new ConsoleAppender(new InfoOnlyLogLayout());
			Logger.getRootLogger().removeAllAppenders();
			Logger.getRootLogger().addAppender(logFileAppender);

			LOGGER.info("Downstreamer" + VERSION);
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
					DownstreamerConverters.convertEqtlToBin(options);
					break;
				case CONVERT_TXT:
					DownstreamerConverters.convertTxtToBin(options);
					break;
				case CONVERT_BIN:
					DownstreamerConverters.convertBinToTxt(options);
					break;
				case CONVERT_GTEX:
					DownstreamerConverters.convertGct(options.getGwasZscoreMatrixPath(), options.getOutputBasePath());
					break;
				case CONVERT_EXP:
					DownstreamerConverters.convertExpressionMatrixToBin(options);
					break;
				case CONVERT_TXT_MERGE:
					DownstreamerConverters.mergeConvertTxt(options);
					break;
				case MERGE_BIN:
					DownstreamerConverters.mergeBinMatrix(options);
					break;
				case FIRST1000:
					First1000qtl.printFirst1000(options);
					break;
				case PTOZSCORE:
					DownstreamerConverters.convertPvalueToZscore(options);
					break;
				case CREATE_EXCEL:
					DownstreamerUtilities.generateExcelFromIntermediates(options);
					break;
				case GET_PATHWAY_LOADINGS:
					DownstreamerUtilities.generatePathwayLoadingExcel(options);
					break;
				case TOP_HITS:
					PruneToIndependentTopHits.prune(options);
					break;
				case STEP1:

					DownstreamerStep1Results step1Res = DownstreamerMainAnalysis.step1(options);
					DownstreamerStep2Results step2Res = null;
					DownstreamerStep3Results step3Res = null;

					if (options.getPathwayDatabases().isEmpty()) {
						LOGGER.info("The analysis will now stop since no pathway databases are provided. Use --mode STEP2 and exactly the same output path and genes file to continue");
						break;
					} else {
						step2Res = DownstreamerMainAnalysis.step2(options, step1Res);
					}

					if (step2Res != null) {
						ExcelWriter writer = null;
						if (options.isSaveOuputAsExcelFiles()) {
							writer = new ExcelWriter(step2Res.getGenePvalues().getColObjects(), options);
							writer.saveStep2Excel(step2Res);
						}
						if (options.isAssignPathwayGenesToCisWindow()) {
							step3Res = DownstreamerMainAnalysis.step3(options);
							if (writer != null) {
								writer.saveStep3Excel(step2Res, step3Res);
							}
						}
					}

					break;

				case STEP2:
					step2Res = DownstreamerMainAnalysis.step2(options, null);

					ExcelWriter writer = null;
					if (options.isSaveOuputAsExcelFiles()) {
						writer = new ExcelWriter(step2Res.getGenePvalues().getColObjects(), options);
						writer.saveStep2Excel(step2Res);
					}
					if (options.isAssignPathwayGenesToCisWindow()) {
						step3Res = DownstreamerMainAnalysis.step3(options);
						if (writer != null) {
							writer.saveStep3Excel(step2Res, step3Res);
						}
					}
					break;
				case CORRELATE_GENES:
					DownstreamerUtilities.correlateGenes(options);
					break;
				case R_2_Z_SCORE:
					DownstreamerUtilities.convertRtoZscore(options);
					break;
				case TRANSPOSE:
					DownstreamerConverters.tranposeBinMatrix(options);
					break;
				case PCA:
					DownstreamerUtilities.doPcaOnBinMatrix(options);
					break;
				case PRIO_GENE_ENRICH:
					TestCoregulationPerformance.testCoreGenePredictionPerformance(options);
					break;
				case INVESTIGATE_NETWORK:
					NetworkProperties.investigateNetwork(options);
					break;
				case GET_NORMALIZED_GENEP:
					DownstreamerUtilities.getNormalizedGwasGenePvalues(options);
					break;
				case REMOVE_CIS_COEXP:
					DownstreamerUtilities.removeLocalGeneCorrelations(options);
					break;
				case SPECIAL:
					CorrelateExpressionToPredictions.run(options);
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

	public static class ThreadErrorHandler implements Thread.UncaughtExceptionHandler {

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

	public static String formatMsForLog(long ms) {
		//return LOG_TIME_FORMAT.format(new Date(ms));
		return DurationFormatUtils.formatDuration(ms, "H:mm:ss.S");
	}

}
