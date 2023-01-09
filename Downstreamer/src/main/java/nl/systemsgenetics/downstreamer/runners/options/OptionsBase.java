package nl.systemsgenetics.downstreamer.runners.options;

import edu.emory.mathcs.utils.ConcurrencyUtils;
import htsjdk.tribble.SimpleFeature;
import org.apache.commons.cli.*;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;

import java.io.File;

public class OptionsBase {
    private static final Logger LOGGER = LogManager.getLogger(OptionsBase.class);
    protected static final Options OPTIONS;

    // Common static variables
    private static int numberOfThreadsToUse = Runtime.getRuntime().availableProcessors();
    private static final SimpleFeature HLA = new SimpleFeature("6", 20000000, 40000000);

    // Common final variables
    private File outputBasePath;
    private File logFile;

    // Common options for all DS modes
    private DownstreamerMode mode;
    private boolean debugMode;
    private File debugFolder;

    static {
        OPTIONS = new Options();

        // mode
        OptionBuilder.withArgName("string");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription(DownstreamerMode.getFullDescriptionString());
        OptionBuilder.withLongOpt("mode");
        OptionBuilder.isRequired();
        OPTIONS.addOption(OptionBuilder.create("m"));

        // debugMode
        OptionBuilder.withArgName("boolean");
        OptionBuilder.withDescription("Activate debug mode. This will result in a more verbose log file and will save many intermediate results to files. Not recommended for large analysis. Not all tools have debug mode.");
        OptionBuilder.withLongOpt("debug");
        OPTIONS.addOption(OptionBuilder.create("d"));

        // outputBasePath
        OptionBuilder.withArgName("path");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("The output path. This is required for logging also when viewing help for a specific tool.");
        OptionBuilder.withLongOpt("output");
        OptionBuilder.isRequired();
        OPTIONS.addOption(OptionBuilder.create("o"));

        // Help
        OptionBuilder.withArgName("boolean");
        OptionBuilder.withDescription("Print help message. Combine with -m for viewing tool specific arguments.");
        OptionBuilder.withLongOpt("help");
        OPTIONS.addOption(OptionBuilder.create("h"));
    }

    public OptionsBase(String[] args) throws ParseException {

        // Parse arguments
        final CommandLineParser parser = new RelaxedParser();
        final CommandLine commandLine = parser.parse(OPTIONS, args, false);

        // Parse the DS mode
        try {
            mode = DownstreamerMode.valueOf(commandLine.getOptionValue('m'));
        } catch (IllegalArgumentException e) {
            LOGGER.error("Invalid mode " + commandLine.getOptionValue('m'));
            LOGGER.error("Please select from the available modes below: ");
            LOGGER.error(DownstreamerMode.getFullDescriptionString());
            throw new ParseException(e.getMessage());
        }

        // Initialize threadpool
        if (commandLine.hasOption('t')) {
            try {
                numberOfThreadsToUse = Integer.parseInt(commandLine.getOptionValue('t'));
                System.setProperty("Djava.util.concurrent.ForkJoinPool.common.parallelism", commandLine.getOptionValue('t'));
                ConcurrencyUtils.setNumberOfThreads(numberOfThreadsToUse);
            } catch (NumberFormatException e) {
                throw new ParseException("Error parsing --threads \"" + commandLine.getOptionValue('t') + "\" is not an int");
            }
        }

        outputBasePath = new File(commandLine.getOptionValue('o'));

        if (outputBasePath.isDirectory()) {
            throw new ParseException("Specified output path is a directory. Please include a prefix for the output files.");
        }

        logFile = new File(outputBasePath + ".log");
        debugFolder = new File(outputBasePath + "_debugFiles");
        debugMode = commandLine.hasOption('d');
    }

	public void printOptions() {

		//TODO subclasses need to overwrite or extend or somthing to print all options relevant
		LOGGER.info("Supplied options:");
		LOGGER.info(" * Mode: " + mode.name());
		LOGGER.info(" * Ouput path: " + outputBasePath.getAbsolutePath());
		LOGGER.info(" * Debug mode: " + (debugMode ? "on (this will result in many intermediate output files)" : "off"));
        LOGGER.info(" * Number of threads to use: " + numberOfThreadsToUse);

	}

    public static void printHelp() {
        HelpFormatter formatter = new HelpFormatter();
        formatter.setWidth(100);
        formatter.printHelp(" ", OPTIONS);
    }

    public static int getNumberOfThreadsToUse() {
        return numberOfThreadsToUse;
    }

    public String getOutputBasePath() {
        return outputBasePath.getPath();
    }

    public File getLogFile() {
        return logFile;
    }

    public DownstreamerMode getMode() {
        return mode;
    }

    public boolean isDebugMode() {
        return debugMode;
    }

    public File getDebugFolder() {
        return debugFolder;
    }
}
