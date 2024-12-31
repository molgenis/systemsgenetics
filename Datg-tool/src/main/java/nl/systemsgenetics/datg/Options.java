package nl.systemsgenetics.datg;

import org.apache.commons.cli.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.File;

public class Options {

    private static final Logger LOGGER = LogManager.getLogger(Options.class);
    protected static final org.apache.commons.cli.Options OPTIONS;

    static {
        OPTIONS = new org.apache.commons.cli.Options();

        // mode
        OptionBuilder.withArgName("string");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription(DatgConvertModes.getFullDescriptionString());
        OptionBuilder.withLongOpt("mode");
        OptionBuilder.isRequired();
        OPTIONS.addOption(OptionBuilder.create("m"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("The output file.");
        OptionBuilder.withLongOpt("output");
        OPTIONS.addOption(OptionBuilder.create("o"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("The input file.");
        OptionBuilder.withLongOpt("input");
        OptionBuilder.isRequired();
        OPTIONS.addOption(OptionBuilder.create("i"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("The type on content on rows IE 'Genes'. This is saved as extra annotation in datg files");
        OptionBuilder.withLongOpt("rowContent");
        OPTIONS.addOption(OptionBuilder.create("rc"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("The type on content on columns IE 'Samples'. This is saved as extra annotation in datg files");
        OptionBuilder.withLongOpt("colContent");
        OPTIONS.addOption(OptionBuilder.create("cc"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("The name of the dataset. This is saved as extra annotation in datg files");
        OptionBuilder.withLongOpt("datasetName");
        OPTIONS.addOption(OptionBuilder.create("dn"));

    }

    private final DatgConvertModes mode;
    private final File outputFile;
    private final File inputFile;
    private final File logFile;
    private final String rowContent;
    private final String colContent;
    private final String datasetName;

    public Options(String[] args) throws Exception {

        final CommandLineParser parser = new PosixParser();
        final CommandLine commandLine = parser.parse(OPTIONS, args, false);

        try {
            String modeString = commandLine.getOptionValue("m").toUpperCase();
            modeString = modeString.replace('-', '_');
            mode = DatgConvertModes.valueOf(modeString);
        } catch (IllegalArgumentException e) {
            throw new ParseException("Error parsing --mode \"" + commandLine.getOptionValue("m") + "\" is not a valid mode");
        }
        String inputArg = commandLine.getOptionValue('i');

        String outputArg = null;
        if(mode != DatgConvertModes.INSPECT){
            outputArg = commandLine.getOptionValue('o');
        }

        switch (mode){
            case DAT_2_DATG:
                if (!inputArg.endsWith(".dat")) {
                    inputArg = inputArg + ".dat";
                }
                //no break intended
            case TXT_2_DATG:
                if(!outputArg.endsWith(".datg")){
                    outputArg = outputArg + ".datg";
                }
                break;
            case DATG_2_TXT:
                if (!inputArg.endsWith(".datg")) {
                    inputArg = inputArg + ".datg";
                }
                if(outputArg.endsWith(".gz")){
                    throw new Exception("Writing as gz file is not supported, it is much faster to just use gzip afterwards");
                }
                if(!outputArg.endsWith(".txt")){
                    outputArg = outputArg + ".txt";
                }
                break;
        }

        if(outputArg == null){
            outputFile = null;
            logFile = null;
        } else {
            outputFile = new File(outputArg);
            logFile = new File((outputArg.endsWith(".datg") ? outputArg.substring(0, outputArg.length() - 5) : outputArg) + ".log");
        }
        inputFile = new File(inputArg);

        rowContent = commandLine.getOptionValue("rowContent","");
        colContent = commandLine.getOptionValue("colContent","");
        datasetName = commandLine.getOptionValue("datasetName","");

    }

    public static void printHelp() {
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp(" ", OPTIONS);
    }

    public void printOptions() {

        LOGGER.info("Supplied options:");
        LOGGER.info(" * Mode: " + mode.name());
        LOGGER.info(" * Input path: " + inputFile.getAbsolutePath());
        if(outputFile != null) {
            LOGGER.info(" * Output path: " + outputFile.getAbsolutePath());
        }
        if(mode != DatgConvertModes.DATG_2_TXT && mode != DatgConvertModes.INSPECT){
            LOGGER.info(" * Row content: " + rowContent);
            LOGGER.info(" * Column content: " + colContent);
            LOGGER.info(" * Dataset name: " + datasetName);
        }


    }

        public DatgConvertModes getMode() {
        return mode;
    }

    public File getOutputFile() {
        return outputFile;
    }

    public File getInputFile() {
        return inputFile;
    }

    public File getLogFile() {
        return logFile;
    }

    public String getDatasetName() {
        return datasetName;
    }

    public String getColContent() {
        return colContent;
    }

    public String getRowContent() {
        return rowContent;
    }
}
