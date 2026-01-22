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

        OptionBuilder.withArgName("string");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("The pattern used to select the files from the input folder. Use brackets to create a capture group that is added to the row names.");
        OptionBuilder.withLongOpt("filePattern");
        OPTIONS.addOption(OptionBuilder.create("fp"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("File with strings per line to grep from row names or single string. Rows with names containing grep string are return by DATG_2_TXT. Case insensitive.");
        OptionBuilder.withLongOpt("rowGrep");
        OPTIONS.addOption(OptionBuilder.create("rg"));



    }

    private final DatgConvertModes mode;
    private final File outputFile;
    private final File inputFile;
    private final File logFile;
    private final String rowContent;
    private final String colContent;
    private final String datasetName;
    private final String filePattern;
    private final File rowGrepFile;
    private final String rowGrepString;
//    private final File colGrepFile;
//    private final String colGrepString;

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
        if(mode != DatgConvertModes.INSPECT && mode != DatgConvertModes.UPGRADE){
            outputArg = commandLine.getOptionValue('o');
        }

        switch (mode){
            case DAT_2_DATG:
                if (!inputArg.endsWith(".dat")) {
                    inputArg = inputArg + ".dat";
                }
                //no break intended
            case ROW_CONCAT:
            case TXT_2_DATG:
                if(!outputArg.endsWith(".datg")){
                    outputArg = outputArg + ".datg";
                }
                break;
            case DATG_2_TXT:
                if(outputArg.endsWith(".gz")){
                    throw new Exception("Writing as gz file is not supported, it is much faster to just use gzip afterwards");
                }
                if(!outputArg.endsWith(".txt")){
                    outputArg = outputArg + ".txt";
                }
                //no break intended
            case INSPECT:
            case UPGRADE:
                if (!inputArg.endsWith(".datg")) {
                    inputArg = inputArg + ".datg";
                }
                break;
        }

        if(outputArg == null){
            outputFile = null;
            if(mode == DatgConvertModes.UPGRADE){
                logFile = new File((inputArg.endsWith(".datg") ? inputArg.substring(0, inputArg.length() - 5) : inputArg) + ".log");
            } else {
                logFile = null;
            }
        } else {
            outputFile = new File(outputArg);
            logFile = new File((outputArg.endsWith(".datg") ? outputArg.substring(0, outputArg.length() - 5) : outputArg) + ".log");
        }
        inputFile = new File(inputArg);

        rowContent = commandLine.getOptionValue("rowContent","");
        colContent = commandLine.getOptionValue("colContent","");
        datasetName = commandLine.getOptionValue("datasetName","");
        filePattern = commandLine.getOptionValue("filePattern",null);

        if(commandLine.hasOption("rowGrep")){
            File file = new File(commandLine.getOptionValue("rowGrep"));
            if(mode != DatgConvertModes.DATG_2_TXT){
                LOGGER.warn("RowGrep is only supported for mode DATG_2_TXT, otherwise it is ignored.");
            }
            if(file.exists()){
                rowGrepFile = file;
                rowGrepString = null;
            } else {
                rowGrepFile = null;
                rowGrepString = commandLine.getOptionValue("rowGrep");
            }
        } else {
            rowGrepFile = null;
            rowGrepString = null;
        }

//        if(commandLine.hasOption("colGrep")){
//            File file = new File(commandLine.getOptionValue("colGrep"));
//            if(mode != DatgConvertModes.DATG_2_TXT){
//                LOGGER.warn("colGrep is only supported for mode DATG_2_TXT, otherwise it is ignored.");
//            }
//            if(file.exists()){
//                colGrepFile = file;
//                colGrepString = null;
//            } else {
//                colGrepFile = null;
//                colGrepString = commandLine.getOptionValue("rowGrep");
//            }
//        } else {
//            colGrepFile = null;
//            colGrepString = null;
//        }




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
        if(mode != DatgConvertModes.DATG_2_TXT && mode != DatgConvertModes.INSPECT && mode != DatgConvertModes.UPGRADE){
            LOGGER.info(" * Row content: " + rowContent);
            LOGGER.info(" * Column content: " + colContent);
            LOGGER.info(" * Dataset name: " + datasetName);
        }
        if(mode != DatgConvertModes.DATG_2_TXT && rowGrepFile != null){
            LOGGER.info(" * Row grep file: " + rowGrepFile.getPath());
        }
        if(mode != DatgConvertModes.DATG_2_TXT && rowGrepString != null){
            LOGGER.info(" * Row grep string: " + rowGrepString);
        }
//        if(mode != DatgConvertModes.DATG_2_TXT && colGrepFile != null){
//            LOGGER.info(" * Column grep file: " + colGrepFile.getPath());
//        }
//        if(mode != DatgConvertModes.DATG_2_TXT && colGrepString != null){
//            LOGGER.info(" * Column grep string: " + colGrepString);
//        }
        if(filePattern != null){
            LOGGER.info(" * File pattern: " + filePattern);
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

    public String getFilePattern() {
        return filePattern;
    }

    public File getRowGrepFile() {
        return rowGrepFile;
    }

    public String getRowGrepString() {
        return rowGrepString;
    }

//    public File getColGrepFile() {
//        return colGrepFile;
//    }
//
//    public String getColGrepString() {
//        return colGrepString;
//    }
}
