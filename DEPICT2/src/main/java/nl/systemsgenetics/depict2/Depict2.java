package nl.systemsgenetics.depict2;

import java.io.File;
import java.io.IOException;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.*;

import nl.systemsgenetics.depict2.development.CorrelateExpressionToPredictions;

import nl.systemsgenetics.depict2.development.First1000qtl;
import nl.systemsgenetics.depict2.io.ExcelWriter;
import nl.systemsgenetics.depict2.runners.*;

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
 * Entry point and error handler for Depict 2
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
	 * Main function for Depict 2 analysis, specific tasks are implemented in
	 * respective runner classes.
	 *
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
			FileAppender logFileAppender = new FileAppender(new SimpleLayout(), options.getLogFile().getCanonicalPath(), options.getMode() != Depict2Mode.STEP1);
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
					Depict2Converters.convertEqtlToBin(options);
					break;
				case CONVERT_TXT:
					Depict2Converters.convertTxtToBin(options);
					break;
				case CONVERT_BIN:
					Depict2Converters.convertBinToTxt(options);
					break;
				case CONVERT_GTEX:
					Depict2Converters.convertGct(options.getGwasZscoreMatrixPath(), options.getOutputBasePath());
					break;
				case CONVERT_EXP:
					Depict2Converters.convertExpressionMatrixToBin(options);
					break;
				case CONVERT_TXT_MERGE:
					Depict2Converters.mergeConvertTxt(options);
					break;
				case MERGE_BIN:
					Depict2Converters.mergeBinMatrix(options);
					break;
				case FIRST1000:
					First1000qtl.printFirst1000(options);
					break;
				case PTOZSCORE:
					Depict2Converters.convertPvalueToZscore(options);
					break;
				case CREATE_EXCEL:
					Depict2Utilities.generateExcelFromIntermediates(options);
					break;
				case GET_PATHWAY_LOADINGS:
					Depict2Utilities.generatePathwayLoadingExcel(options);
					break;
				case TOP_HITS:
					PruneToIndependentTopHits.prune(options);
					break;
				case STEP1:

					Depict2Step1Results step1Res = Depict2MainAnalysis.step1(options);
					Depict2Step2Results step2Res = null;
					Depict2Step3Results step3Res = null;

					if (options.getPathwayDatabases().isEmpty()) {
						LOGGER.info("The analysis will now stop since no pathway databases are provided. Use --mode STEP2 and exactly the same output path and genes file to continue");
						break;
					} else {
						step2Res = Depict2MainAnalysis.step2(options, step1Res);
					}

					if (step2Res != null) {
						ExcelWriter writer = null;
						if (options.isSaveOuputAsExcelFiles()) {
							writer = new ExcelWriter(step2Res.getGenePvalues().getColObjects(), options);
							writer.saveStep2Excel(step2Res);
						}
						if (options.isAssignPathwayGenesToCisWindow()) {
							step3Res = Depict2MainAnalysis.step3(options);
							if (writer != null) {
								writer.saveStep3Excel(step2Res, step3Res);
							}
						}
					}

					break;

				case STEP2:
					step2Res = Depict2MainAnalysis.step2(options, null);

					ExcelWriter writer = null;
					if (options.isSaveOuputAsExcelFiles()) {
						writer = new ExcelWriter(step2Res.getGenePvalues().getColObjects(), options);
						writer.saveStep2Excel(step2Res);
					}
					if (options.isAssignPathwayGenesToCisWindow()) {
						step3Res = Depict2MainAnalysis.step3(options);
						if (writer != null) {
							writer.saveStep3Excel(step2Res, step3Res);
						}
					}
					break;
				case CORRELATE_GENES:
					Depict2Utilities.correlateGenes(options);
					break;
				case R_2_Z_SCORE:
					Depict2Utilities.convertRtoZscore(options);
					break;
				case TRANSPOSE:
					Depict2Converters.tranposeBinMatrix(options);
					break;
				case PCA:
					Depict2Utilities.doPcaOnBinMatrix(options);
					break;
				case CORE_GENE_AUC:
					TestCoregulationPerformance.testCoreGenePredictionPerformance(options);
					break;
				case INVESTIGATE_NETWORK:
					NetworkProperties.investigateNetwork(options);
					break;
				case GET_NORMALIZED_GENEP:
					Depict2Utilities.getNormalizedGwasGenePvalues(options);
					break;
				case REMOVE_LOCAL_COR:
					Depict2Utilities.removeLocalGeneCorrelations(options);
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
