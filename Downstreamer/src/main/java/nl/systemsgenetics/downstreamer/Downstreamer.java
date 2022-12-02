package nl.systemsgenetics.downstreamer;

import nl.systemsgenetics.downstreamer.io.ExcelWriter;
import nl.systemsgenetics.downstreamer.pathway.PredictedPathwayAnnotations;
import nl.systemsgenetics.downstreamer.runners.options.*;
import nl.systemsgenetics.downstreamer.runners.*;
import org.apache.commons.cli.ParseException;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.lang.time.DurationFormatUtils;
import org.apache.log4j.*;

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

    private static final Logger LOGGER = Logger.getLogger(Downstreamer.class);
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

        OptionsBase options;

        // If no arguments are provided print help
        if (args.length == 0) {
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
            // Append to the existing logfile
            FileAppender logFileAppender = new FileAppender(new SimpleLayout(), options.getLogFile().getCanonicalPath(), true);
            ConsoleAppender logConsoleInfoAppender = new ConsoleAppender(new InfoOnlyLogLayout());
            Logger.getRootLogger().removeAllAppenders();
            Logger.getRootLogger().addAppender(logFileAppender);

            // Create a seperator in the log, so seperate runs are easier to spot
            LOGGER.info("\n");
            LOGGER.info(StringUtils.repeat("-", 80));
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
                case REGRESS: DownstreamerRegressionEnginePColt.run(new OptionsModeRegress(args)); break;

                // Converter utils that share OptionsModeConvert()
                case CONVERT_TXT: DownstreamerConverters.convertTxtToBin(new OptionsModeConvert(args)); break;
                case CONVERT_BIN: DownstreamerConverters.convertBinToTxt(new OptionsModeConvert(args)); break;
                case CONVERT_MERGE_BIN: DownstreamerConverters.mergeBinMatrix(new OptionsModeConvert(args)); break;
                case CONVERT_EXP: DownstreamerConverters.convertExpressionMatrixToBin(new OptionsModeConvert(args)); break;
                case CONVERT_DAT_TO_DATG: DownstreamerConverters.convertDatToDatg(new OptionsModeConvert(args)); break;
                case CONVERT_GTEX: DownstreamerConverters.convertGct(new OptionsModeConvert(args)); break;
                case CONVERT_EQTL: DownstreamerConverters.convertEqtlToBin(new OptionsModeConvert(args)); break;
                case CONVERT_PTOZSCORE: DownstreamerConverters.convertPvalueToZscore(new OptionsModeConvert(args)); break;
                case PREPARE_GWAS: DownstreamerConverters.prepareGwasSummaryStatistics(new OptionsModeConvert(args)); break;
                case PREPARE_GWAS_MERGE: DownstreamerConverters.prepareGwasSummaryStatisticsMerge(new OptionsModeConvert(args)); break;
                case MATRIX_SUBSET: DownstreamerConverters.subsetMatrix(new OptionsModeConvert(args)); break;
                case MATRIX_TRANSPOSE: DownstreamerConverters.tranposeBinMatrix(new OptionsModeConvert(args));break;

                // Utilities related to making the co-regulation network. Share OptionsModeCoreg()
                case COREG_CORRELATE_GENES: CoregulationUtilities.coregCorrelateGenes(new OptionsModeCoreg(args)); break;
                case COREG_RTOZSCORE: CoregulationUtilities.coregConvertRtoZscore(new OptionsModeCoreg(args)); break;
                case COREG_PCA: CoregulationUtilities.coregDoPcaOnBinMatrix(new OptionsModeCoreg(args)); break;
                case COREG_REMOVE_CIS_COEXP: CoregulationUtilities.coregRemoveLocalGeneCorrelations(new OptionsModeCoreg(args));break;
                case COREG_INVESTIGATE_NETWORK: CoregulationUtilities.coregInvestigateNetwork(new OptionsModeCoreg(args)); break;

                // Legacy modes, These will be either replaced or substantially updated
                case CREATE_EXCEL:
                    DownstreamerUtilities.generateExcelFromIntermediates(new DownstreamerOptionsDeprecated(args));
                    break;
                case PRIO_GENE_ENRICH:
                    PathwayDatabaseEnrichments.testPredictionPerformance(new DownstreamerOptionsDeprecated(args));
                    break;
                case EXPAND_PATHWAYS:
                    PredictedPathwayAnnotations.expandAnnotations(new DownstreamerOptionsDeprecated(args));
                    break;
                case STEP1:
                    optionsDeprecated = new DownstreamerOptionsDeprecated(args);
                    final DownstreamerStep1Results step1Res = DownstreamerMainAnalysis.step1(optionsDeprecated);
                    final DownstreamerStep2Results step2Res;

                    if (optionsDeprecated.getPathwayDatabases().isEmpty()) {
                        LOGGER.info("The analysis will now stop since no pathway databases are provided. Use --mode STEP2 and exactly the same output path and genes file to continue");
                        break;
                    } else {
                        step2Res = DownstreamerMainAnalysis.step2(optionsDeprecated, step1Res);
                    }

                    if (step2Res != null) {
                        ExcelWriter writer = null;
                        if (optionsDeprecated.isSaveOuputAsExcelFiles()) {
                            writer = new ExcelWriter(step2Res.getGenePvalues().getColObjects(), optionsDeprecated);
                            writer.saveStep2Excel(step2Res);
                        }

                        if (optionsDeprecated.isAssignPathwayGenesToCisWindow()) {
                            final DownstreamerStep3Results step3Res = DownstreamerMainAnalysis.step3(optionsDeprecated);
                            if (writer != null) {
                                writer.saveStep3Excel(step2Res, step3Res);
                            }
                        }
                    }

                    break;

                case STEP2:
                    optionsDeprecated = new DownstreamerOptionsDeprecated(args);
                    step2Res = DownstreamerMainAnalysis.step2(optionsDeprecated, null);

                    ExcelWriter writer = null;
                    if (optionsDeprecated.isSaveOuputAsExcelFiles()) {
                        writer = new ExcelWriter(step2Res.getGenePvalues().getColObjects(), optionsDeprecated);
                        writer.saveStep2Excel(step2Res);
                    }
                    if (optionsDeprecated.isAssignPathwayGenesToCisWindow()) {
                        final DownstreamerStep3Results step3Res = DownstreamerMainAnalysis.step3(optionsDeprecated);
                        if (writer != null) {
                            writer.saveStep3Excel(step2Res, step3Res);
                        }
                    }
                    break;
                default:
                    throw new Exception("Mode not yet converted or is deprecated.");
            }
        } catch (ParseException e) {
            System.err.println("Problem running mode: " + options.getMode());
            System.err.println("Error message: " + e.getMessage());
            LOGGER.fatal("Error: " + e.getMessage(), e);
            helpProvider.printHelp();
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
        LOGGER.info(StringUtils.repeat("-", 80));
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
