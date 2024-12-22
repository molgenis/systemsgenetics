package nl.systemsgenetics.datg;

import com.google.common.io.Files;
import com.opencsv.*;
import org.apache.commons.cli.ParseException;
import org.apache.commons.math3.util.Precision;
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
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRow;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRowCompressedReader;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRowCompressedWriter;

import java.io.*;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.ResourceBundle;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

public class DatgConverter {

    public static final String VERSION = ResourceBundle.getBundle("verion").getString("application.version");
    public static final DateFormat DATE_TIME_FORMAT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
    public static final DateFormat TIME_FORMAT = new SimpleDateFormat("HH:mm:ss");

    private static final Logger LOGGER = LogManager.getLogger(DatgConverter.class);
    private static final String HEADER
            = "  /---------------------------------------\\\n"
            + "  |            Datg Converter             |\n"
            + "  |                                       |\n"
            + "  |  University Medical Center Groningen  |\n"
            + "  \\---------------------------------------/";


    public static void main(String[] args) throws InterruptedException {
        System.out.println();
        System.out.println(HEADER);
        System.out.println();
        System.out.println("       --- Version: " + VERSION + " ---");
        System.out.println();
        System.out.println("More information: https://github.com/molgenis/systemsgenetics/wiki/");
        System.out.println();

        Date currentDataTime = new Date();
        String startDateTime = DATE_TIME_FORMAT.format(currentDataTime);
        System.out.println("Current date and time: " + startDateTime);
        System.out.println();

        System.out.flush(); //flush to make sure header is before errors
        Thread.sleep(25); //Allows flush to complete

        final Options options;
        // Parse general options
        try {
            options = new Options(args);
        } catch (ParseException ex) {
            System.err.println("Error parsing commandline: " + ex.getMessage());
            Options.printHelp();
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
            initializeLoggers(loggingLevel, options.getLogFile(), startDateTime);

        } catch (IOException e) {
            System.err.println("Failed to create logger: " + e.getMessage());
            System.exit(1);
        }

        options.printOptions();

        try{
            switch (options.getMode()){

                case TXT_2_DATG:
                    //DoubleMatrixDataset.loadDoubleTextData(options.getInputFile().getAbsolutePath());

                    convertTxt2Datg(options);

                    break;
                case DATG_2_TXT:
                    DoubleMatrixDatasetRowCompressedReader datgData = new DoubleMatrixDatasetRowCompressedReader(options.getInputFile().getAbsolutePath());

                    LOGGER.info("Input .datg file contains " + datgData.getNumberOfColumns() + " columns and " + datgData.getNumberOfRows());
                    LOGGER.info("The dataset is named: '" + datgData.getDatasetName() + "' and was created on: " +  DATE_TIME_FORMAT.format(datgData.getTimestampCreated()));

                    final CSVWriter outputWriter;
                    outputWriter = new CSVWriter(new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(options.getOutputFile())))), '\t', '\0', '\0', "\n");

                    final String[] outputLine = new String[datgData.getNumberOfColumns() + 1];
                    int c = 0;
                    outputLine[c++] = "";
                    for(String column : datgData.getColumnIdentifiers()){
                        outputLine[c++] = column;
                    }

                    outputWriter.writeNext(outputLine);

                    for(DoubleMatrixDatasetRow rowData : datgData){

                        c = 0;
                        outputLine[c++] = rowData.getRowName();
                        double[] rowContent = rowData.getData();
                        for(int i = 0 ; i < datgData.getNumberOfColumns() ; ++i){
                            outputLine[c++] = String.valueOf(rowContent[i]);
                        }

                        outputWriter.writeNext(outputLine);


                    }
                    outputWriter.close();

                    break;
                case DAT_2_DATG:

                    //reading the whole dataset is a bit inefficient but this is only a legacy function convert old files

                    DoubleMatrixDataset<String, String> data = DoubleMatrixDataset.loadDoubleBinaryData(options.getInputFile().getAbsolutePath());
                    LOGGER.info("Input .dat file contains " + data.columns() + " columns and " + data.rows());
                    data.saveBinary(options.getOutputFile().getAbsolutePath(), options.getDatasetName(), options.getRowContent(), options.getColContent());

                    break;
            }
        } catch (IOException e) {
                System.err.println("Problem running mode: " + options.getMode());
                LOGGER.fatal("Error: ", e);
                System.exit(1);
            }

        currentDataTime = new Date();
        LOGGER.info("Conversion completed at: " + DATE_TIME_FORMAT.format(currentDataTime));

    }

    private static void convertTxt2Datg(Options options) throws IOException {
        final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
        final CSVReader reader = new CSVReaderBuilder((new BufferedReader(new InputStreamReader((new GZIPInputStream(new FileInputStream(options.getInputFile()))))))).withCSVParser(parser).build();

        //First nextLine contains the header
        String[] nextLine = reader.readNext();

        ArrayList<String> colnames = new ArrayList<>(nextLine.length-1);
        for(int i = 1 ; i < nextLine.length ; i++){
            colnames.add(nextLine[i]);
        }

        DoubleMatrixDatasetRowCompressedWriter datgWriter = new DoubleMatrixDatasetRowCompressedWriter(options.getOutputFile().getAbsolutePath(), colnames, options.getDatasetName(), options.getRowContent(), options.getColContent());

        double[] rowData = new double[colnames.size()];
        while ((nextLine = reader.readNext()) != null) {

            for(int i = 1 ; i < nextLine.length ; i++){
                rowData[i] = Double.parseDouble(nextLine[i]);
            }

            datgWriter.addRow(nextLine[0], rowData);
        }

        datgWriter.close();
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
        LOGGER.info("Datg converter " + VERSION);
        LOGGER.info("Current date and time: " + startDateTime);

        rootLoggerConfig.addAppender(stdOut, Level.INFO, null);
    }


}
