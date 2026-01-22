package nl.systemsgenetics.datg;

import com.opencsv.*;
import net.jpountz.lz4.LZ4FrameOutputStream;
import org.apache.commons.cli.ParseException;
import org.apache.commons.collections.CollectionUtils;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRow;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRowCompressedReader;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRowCompressedWriter;
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
import org.ahocorasick.trie.*;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

public class DatgConverter {

    public static final String VERSION = ResourceBundle.getBundle("verion").getString("application.version");
    public static final DateFormat DATE_TIME_FORMAT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

    //This is set to true by testng to throw error instead of catch it for nice user output.
    protected static boolean TESTNG_MODE = false;

    private static final Logger LOGGER = LogManager.getLogger(DatgConverter.class);
    private static final String HEADER
            = "  /----------------------------------------\\\n"
            + "  |             Datg Converter             |\n"
            + "  |                                        |\n"
            + "  | University Medical Center Groningen \uD83D\uDC38 |\n"
            + "  \\----------------------------------------/";


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
        } catch (Exception ex) {
            System.err.println("Error parsing commandline: " + ex.getMessage());
            return;
        }

        if (options.getOutputFile() != null) {
            // If the output folder does not exist, create it
            if (options.getLogFile().getParentFile() != null && !options.getLogFile().getParentFile().isDirectory()) {
                if (!options.getLogFile().getParentFile().mkdirs()) {
                    System.err.println("Failed to create output folder: " + options.getLogFile().getParent());
                    System.exit(1);
                }
            }
        }

        // Initialize logfile
        try {

            Level loggingLevel = Level.INFO;
            initializeLoggers(loggingLevel, options.getLogFile(), startDateTime);

        } catch (Exception e) {
            System.err.println("Failed to create logger: " + e.getMessage());
            if (TESTNG_MODE) {
                throw new RuntimeException(e);
            } else {
                System.exit(1);
            }
        }

        options.printOptions();

        try {
            switch (options.getMode()) {

                case TXT_2_DATG:
                    //DoubleMatrixDataset.loadDoubleTextData(options.getInputFile().getAbsolutePath());

                    convertTxt2Datg(options);

                    break;
                case DATG_2_TXT:
                    convertDatg2Txt(options);

                    break;
                case DAT_2_DATG:

                    //reading the whole dataset is a bit inefficient but this is only a legacy function to convert old files
                {
                    DoubleMatrixDataset<String, String> data = DoubleMatrixDataset.loadDoubleBinaryData(options.getInputFile().getAbsolutePath());
                    LOGGER.info("Input .dat file contains " + data.columns() + " columns and " + data.rows() + " rows");
                    data.saveBinary(options.getOutputFile().getAbsolutePath(), options.getDatasetName(), options.getRowContent(), options.getColContent());
                }
                break;
                case INSPECT:

                    inspect(options);
                    break;

                case TEST:

                    testFunction(options);
                    break;

                case ROW_CONCAT:

                    rowConcat(options);
                    break;

                case UPGRADE:

                    upgrade(options);
                    break;

            }
        } catch (Exception e) {
            System.err.println("Problem running mode: " + options.getMode());
            LOGGER.fatal("Error: ", e);
            if (TESTNG_MODE) {
                throw new RuntimeException(e);
            } else {
                System.exit(1);
            }
        }

        currentDataTime = new Date();
        LOGGER.info("Conversion completed at: " + DATE_TIME_FORMAT.format(currentDataTime));

    }

    private static void upgrade(Options options) throws IOException {

        String inputMatrix = options.getInputFile().getAbsolutePath();

        if (inputMatrix.endsWith(".datg")) {
            inputMatrix = inputMatrix.substring(0, inputMatrix.length() - 5);
        }

        File originalDatg = new File(inputMatrix + ".datg");
        File originalRow = new File(inputMatrix + ".rows.txt.gz");
        File originalCol = new File(inputMatrix + ".cols.txt.gz");

        File workdir = originalDatg.getParentFile();
        File tmpDir = new File(workdir, "tmpdir_" + String.valueOf(Math.abs(new Random().nextInt())));
        if (!tmpDir.mkdirs()) {
            throw new IOException("Failed to create tmp dir: " + tmpDir);
        }

        LOGGER.info("Tmp dir: " + tmpDir.getAbsolutePath());

        File originalDatTmp = new File(tmpDir, originalDatg.getName());
        File originalRowTmp = new File(tmpDir, originalRow.getName());
        File originalColTmp = new File(tmpDir, originalCol.getName());

        Files.move(originalDatg.toPath(), originalDatTmp.toPath());
        Files.move(originalRow.toPath(), originalRowTmp.toPath());
        Files.move(originalCol.toPath(), originalColTmp.toPath());

        DoubleMatrixDatasetRowCompressedReader originalReader = new DoubleMatrixDatasetRowCompressedReader(originalDatTmp.getAbsolutePath());

        DoubleMatrixDatasetRowCompressedWriter upgradedWriter = new DoubleMatrixDatasetRowCompressedWriter(inputMatrix, originalReader.getColumnIdentifiers(), originalReader.getDatasetName(), originalReader.getDataOnRows(), originalReader.getDataOnCols());

        for (DoubleMatrixDatasetRow row : originalReader) {
            upgradedWriter.addRow(row.getRowName(), row.getData());
        }

        upgradedWriter.close();


        DoubleMatrixDatasetRowCompressedReader updatedReader = new DoubleMatrixDatasetRowCompressedReader(inputMatrix);


        try {
           validate(originalReader, updatedReader);
        } catch (Exception ex) {
            LOGGER.error("Error upgrading the matrix: {}", ex.getMessage());
            new File(inputMatrix).delete();
            LOGGER.error("Please use backup of original data: {}", originalDatTmp.getAbsolutePath());
            return;
        }

        //Here we only get if the comparison is succesfull
        LOGGER.info("New file is identical to original, removing tmp files");

        originalDatTmp.delete();
        originalRowTmp.delete();
        originalColTmp.delete();

        if (tmpDir.listFiles().length == 0) {
            tmpDir.delete();
        }

    }

    private static void validate(DoubleMatrixDatasetRowCompressedReader originalReader, DoubleMatrixDatasetRowCompressedReader updatedReader) throws Exception {

        if (updatedReader.getNumberOfRows() != originalReader.getNumberOfRows()) {
            throw new Exception("Number rows not equal");
        }

        if (updatedReader.getNumberOfColumns() != originalReader.getNumberOfColumns()) {
            throw new Exception("Number columns not equal");
        }

        if (!originalReader.getDatasetName().equals(updatedReader.getDatasetName()) ||
                !originalReader.getDataOnRows().equals(updatedReader.getDataOnRows()) ||
                !originalReader.getDataOnCols().equals(updatedReader.getDataOnCols())) {
            throw new Exception("Meta data not equal");
        }

        Iterator<String> originalColumniterator = originalReader.getColumnIdentifiers().iterator();
        Iterator<String> updatedColumniterator = updatedReader.getColumnIdentifiers().iterator();

        while (originalColumniterator.hasNext()) {
            if(!originalColumniterator.next().equals(updatedColumniterator.next())) {
                throw new Exception("Different column names");
            }
        }

        for (DoubleMatrixDatasetRow row : originalReader) {
            double[] oldData = row.getData();
            double[] newData = updatedReader.loadSingleRow(row.getRowName());
            for(int i = 0 ; i < oldData.length ; i++) {
                if(oldData[i] != newData[i]) {
                    throw new Exception("Data difference");
                }
            }
        }

    }

    private static void rowConcat(Options options) throws IOException {

        //Input is folder with files to parse
        final File inputFolder = options.getInputFile();
        final Pattern filePattern = Pattern.compile(options.getFilePattern());

        if (!inputFolder.exists()) {
            throw new IOException("Input folder not found at: " + inputFolder.getPath());
        }

        if (!inputFolder.isDirectory()) {
            throw new IOException("Input folder is not a directory at: " + inputFolder.getPath());
        }

        final File[] files = inputFolder.listFiles();

        if (files == null || files.length == 0) {
            throw new IOException("Input folder is empty at: " + inputFolder.getPath());
        }

        final boolean addToRowNames;
        if (filePattern.matcher("").groupCount() == 1) {
            LOGGER.info("Found capturing group in file pattern, will be added to row names.");
            addToRowNames = true;
        } else {
            LOGGER.info("No capturing group in file pattern, rows will keep original names.");
            addToRowNames = false;
        }

        final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();

        ArrayList<String> colnames = null;
        DoubleMatrixDatasetRowCompressedWriter datgWriter = null;

        int includedFilesCount = 0;
        for (File inputFile : files) {
            final Matcher fileMatcher = filePattern.matcher(inputFile.getName());
            if (!fileMatcher.find()) {
                continue;
            }
            final String rowNameAppend = addToRowNames ? fileMatcher.group(1) : null;

            LOGGER.info(String.format("Adding file %s", inputFile.getName()));

            final CSVReader reader;
            if (inputFile.getName().endsWith(".gz")) {
                reader = new CSVReaderBuilder((new BufferedReader(new InputStreamReader((new GZIPInputStream(new FileInputStream(inputFile))))))).withCSVParser(parser).build();
            } else {
                reader = new CSVReaderBuilder((new BufferedReader(new InputStreamReader((new FileInputStream(inputFile)))))).withCSVParser(parser).build();
            }

            //First nextLine contains the header
            String[] nextLine = reader.readNext();

            if (nextLine == null) {
                throw new IOException("File is empty");
            }

            if (includedFilesCount == 0) {
                if (nextLine.length < 2) {
                    throw new IOException("No columns found in: " + inputFile);
                }
                colnames = new ArrayList<>(nextLine.length - 1);
                for (int i = 1; i < nextLine.length; i++) {
                    colnames.add(nextLine[i]);
                }
            } else {
                if (nextLine.length - 1 != colnames.size()) {
                    throw new IOException("Number of columns not identical to first file.");
                }
                for (int i = 1; i < nextLine.length; i++) {
                    if (!colnames.get(i - 1).equals(nextLine[i])) {
                        throw new IOException("All files need identical columns in same order.");
                    }
                }
            }

            if (includedFilesCount == 0) {

                datgWriter = new DoubleMatrixDatasetRowCompressedWriter(
                        options.getOutputFile().getAbsolutePath(),
                        colnames,
                        options.getDatasetName(),
                        options.getRowContent(),
                        options.getColContent()
                );
            }

            double[] rowData = new double[colnames.size()];
            while ((nextLine = reader.readNext()) != null) {

                String rowName = addToRowNames ? nextLine[0] + "-" + rowNameAppend : nextLine[0];

                for (int i = 1; i < nextLine.length; i++) {

                    //Some manual extra to parse double to be compatible with other tools
                    if (nextLine[i].isEmpty()) {
                        rowData[i - 1] = Double.NaN;
                    } else if (nextLine[i].equals("inf")) {
                        rowData[i - 1] = Double.POSITIVE_INFINITY;
                    } else if (nextLine[i].equals("-inf")) {
                        rowData[i - 1] = Double.NEGATIVE_INFINITY;
                    } else {
                        rowData[i - 1] = Double.parseDouble(nextLine[i]);
                    }

                }

                datgWriter.addRow(rowName, rowData);
            }


            includedFilesCount++;
        }

        LOGGER.info(String.format("Included %d files", includedFilesCount));

        if (datgWriter != null) {
            datgWriter.close();
        }

    }

    private static void inspect(Options options) throws IOException {
        DoubleMatrixDatasetRowCompressedReader datgData = new DoubleMatrixDatasetRowCompressedReader(options.getInputFile().getAbsolutePath());

        LOGGER.info(String.format("The dataset is named: '%s'", datgData.getDatasetName()));
        LOGGER.info(String.format("The dataset was created on: %s", DATE_TIME_FORMAT.format(new Date(datgData.getTimestampCreated() * 1000))));
        LOGGER.info(String.format("Number of rows: %d content on rows: '%s'", datgData.getNumberOfRows(), datgData.getDataOnRows()));
        LOGGER.info(String.format("Number of columns: %d content on columns: '%s'", datgData.getNumberOfColumns(), datgData.getDataOnCols()));
        LOGGER.info(String.format("Compression algorithm: '%s'", datgData.isLz4FrameCompression() ? "LZ4 frames" : "jpountz.lz4 block variant"));

        LOGGER.info("");
        LOGGER.info("Top left corner:");
        LOGGER.info("");

        final int colsToShow = Math.min(datgData.getNumberOfColumns(), 4);
        final int rowsToShow = Math.min(datgData.getNumberOfRows(), 10);

        Iterator<String> columnsItt = datgData.getColumnIdentifiers().iterator();
        Iterator<String> rowsItt = datgData.getRowIdentifiers().iterator();

        StringBuilder printRow = new StringBuilder(String.format("%21s  ", ""));
        for (int c = 0; c < colsToShow; c++) {

            String colName = columnsItt.next();
            if (colName.length() > 20) {
                colName = colName.substring(0, 8) + "..." + colName.substring(colName.length() - 8);
            }

            printRow.append(String.format("%21s │", colName));
        }

        LOGGER.info(printRow.toString());

        final DecimalFormat normalFormatter = new DecimalFormat("0.####");
        final DecimalFormat scientificFormatter = new DecimalFormat("0.####E0");

        final String rowSepBlock = "══════════════════════";
        printRow = new StringBuilder("                      ╔");
        for (int c = 0; c < colsToShow; c++) {
            printRow.append(rowSepBlock);
            if (c == colsToShow - 1) {
                printRow.append('╡');
            } else {
                printRow.append('╪');
            }
        }
        LOGGER.info(printRow.toString());

        for (int r = 0; r < rowsToShow; r++) {
            String rowName = rowsItt.next();
            if (rowName.length() > 20) {
                rowName = rowName.substring(0, 8) + "..." + rowName.substring(rowName.length() - 8);
            }
            printRow = new StringBuilder(String.format("%21s ║", rowName));
            double[] rowData = datgData.loadSingleRow(r);
            for (int c = 0; c < colsToShow; c++) {
                double value = rowData[c];
                double valueAbs = Math.abs(value);
                String valueFormatted;
                if (value == 0 || (valueAbs >= 0.001 && valueAbs < 10000)) {
                    valueFormatted = normalFormatter.format(value);
                } else {
                    valueFormatted = scientificFormatter.format(value);
                }
                printRow.append(String.format("%21s │", valueFormatted));
            }
            LOGGER.info(printRow.toString());
        }

        LOGGER.info("");
        LOGGER.info("Note values are rounded and long names shortened with ... in the middle.");
        LOGGER.info("");
    }

    private static void convertDatg2Txt(Options options) throws IOException {
        DoubleMatrixDatasetRowCompressedReader datgData = new DoubleMatrixDatasetRowCompressedReader(options.getInputFile().getAbsolutePath());

        LOGGER.info("Input .datg file contains " + datgData.getNumberOfColumns() + " columns " + (datgData.getDataOnCols().isEmpty() ? "" : "(" + datgData.getDataOnCols() + ") ") + "and " + datgData.getNumberOfRows() + " rows" + (datgData.getDataOnRows().isEmpty() ? "" : " (" + datgData.getDataOnRows() + ")"));
        LOGGER.info("The dataset is named: '" + datgData.getDatasetName() + "' and was created on: " + DATE_TIME_FORMAT.format(new Date(datgData.getTimestampCreated() * 1000)));

        final CSVWriter outputWriter;
        outputWriter = new CSVWriter(new OutputStreamWriter(new BufferedOutputStream(new FileOutputStream(options.getOutputFile())), StandardCharsets.UTF_8), '\t', '\0', '\0', "\n");

        final String[] outputLine = new String[datgData.getNumberOfColumns() + 1];
        int c = 0;
        outputLine[c++] = "";
        for (String column : datgData.getColumnIdentifiers()) {
            outputLine[c++] = column;
        }

        outputWriter.writeNext(outputLine);


        if (options.getRowGrepFile() != null || options.getRowGrepString() != null) {

            final Trie.TrieBuilder grepTrieBuilder = Trie.builder().ignoreCase();
            int numberOfGrepTerms = 0;

            if (options.getRowGrepFile() != null) {

                if (!options.getRowGrepFile().exists()) {
                    throw new IOException("Row grep file does not exist");
                }

                final BufferedReader reader;
                if (options.getRowGrepFile().getName().endsWith(".gz")) {
                    reader = (new BufferedReader(new InputStreamReader((new GZIPInputStream(new FileInputStream(options.getRowGrepFile()))))));
                } else {
                    reader = (new BufferedReader(new InputStreamReader((new FileInputStream(options.getRowGrepFile())))));
                }


                String line;
                while ((line = reader.readLine()) != null) {
                    grepTrieBuilder.addKeyword(line);
                    numberOfGrepTerms++;
                }


            } else {

                //Using tree for simply query is a bit overkill.
                grepTrieBuilder.addKeyword(options.getRowGrepString());
                numberOfGrepTerms++;

            }

            LOGGER.info("Number of row grep terms: " + numberOfGrepTerms);
            final Trie rowGrepTrie = grepTrieBuilder.build();

            int countRowMatched = 0;
            for (Map.Entry<String, Integer> rowEntry : datgData.getRowMap().entrySet()) {
                String rowName = rowEntry.getKey();
                if (rowGrepTrie.containsMatch(rowName)) {
                    ++countRowMatched;

                    double[] rowContent = datgData.loadSingleRow(rowEntry.getValue());

                    c = 0;
                    outputLine[c++] = rowName;
                    for (int i = 0; i < datgData.getNumberOfColumns(); ++i) {
                        outputLine[c++] = String.valueOf(rowContent[i]);
                    }

                    outputWriter.writeNext(outputLine);

                }
            }
            LOGGER.info("Number of rows matching to grep: " + countRowMatched);

        } else {
            for (DoubleMatrixDatasetRow rowData : datgData) {

                c = 0;
                outputLine[c++] = rowData.getRowName();
                double[] rowContent = rowData.getData();
                for (int i = 0; i < datgData.getNumberOfColumns(); ++i) {
                    outputLine[c++] = String.valueOf(rowContent[i]);
                }

                outputWriter.writeNext(outputLine);

            }

        }

        outputWriter.close();

    }

    private static void testFunction(Options options) throws IOException {

        Random random = new Random(42);
        DoubleMatrixDataset<String, String> testData = new DoubleMatrixDataset<>(20000, 50000);
        for (int r = 0; r < testData.rows(); ++r) {
            for (int c = 0; c < testData.columns(); ++c) {
                testData.setElementQuick(r, c, random.nextDouble());
            }
        }

        LOGGER.info("T1 " + DATE_TIME_FORMAT.format(new Date()));

        DoubleMatrixDatasetRowCompressedWriter.saveDataset(options.getOutputFile().getPath(), testData);

        LOGGER.info("T2 " + DATE_TIME_FORMAT.format(new Date()));

        DoubleMatrixDataset<String, String> tmp = new DoubleMatrixDatasetRowCompressedReader(options.getOutputFile().getPath()).loadFullDataset();

        LOGGER.info("T3 " + DATE_TIME_FORMAT.format(new Date()));

        testData.saveBinaryOldFormat(options.getOutputFile().getPath() + "old");

        LOGGER.info("T4 " + DATE_TIME_FORMAT.format(new Date()));

        DoubleMatrixDataset.loadDoubleBinaryData(options.getOutputFile().getPath() + "old");

        LOGGER.info("T5 " + DATE_TIME_FORMAT.format(new Date()));

    }

    private static void convertTxt2Datg(Options options) throws IOException {

        if (!options.getInputFile().exists()) {
            throw new IOException("Input file not found at: " + options.getInputFile());
        }

        final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
        final CSVReader reader;
        if (options.getInputFile().getName().endsWith(".gz")) {
            reader = new CSVReaderBuilder((new BufferedReader(new InputStreamReader((new GZIPInputStream(new FileInputStream(options.getInputFile()))))))).withCSVParser(parser).build();
        } else {
            reader = new CSVReaderBuilder((new BufferedReader(new InputStreamReader((new FileInputStream(options.getInputFile())))))).withCSVParser(parser).build();
        }


        //First nextLine contains the header
        String[] nextLine = reader.readNext();

        ArrayList<String> colnames = new ArrayList<>(nextLine.length - 1);
        for (int i = 1; i < nextLine.length; i++) {
            colnames.add(nextLine[i]);
        }

        if (colnames.isEmpty()) {
            throw new IOException("No columns found in: " + options.getInputFile());
        }


        DoubleMatrixDatasetRowCompressedWriter datgWriter =
                new DoubleMatrixDatasetRowCompressedWriter(
                        options.getOutputFile().getAbsolutePath(),
                        colnames,
                        options.getDatasetName(),
                        options.getRowContent(),
                        options.getColContent()
                );

        double[] rowData = new double[colnames.size()];
        while ((nextLine = reader.readNext()) != null) {

            for (int i = 1; i < nextLine.length; i++) {

                //Some manual extra to parse double to be compatible with other tools
                if (nextLine[i].isEmpty()) {
                    rowData[i - 1] = Double.NaN;
                } else if (nextLine[i].equals("inf")) {
                    rowData[i - 1] = Double.POSITIVE_INFINITY;
                } else if (nextLine[i].equals("-inf")) {
                    rowData[i - 1] = Double.NEGATIVE_INFINITY;
                } else {
                    rowData[i - 1] = Double.parseDouble(nextLine[i]);
                }

            }

            datgWriter.addRow(nextLine[0], rowData);
        }

        datgWriter.close();
    }

    public static void initializeLoggers(Level loggingLevel, File logFile, String startDateTime) throws IOException {
        Configurator.setRootLevel(loggingLevel);
        LoggerContext context = LoggerContext.getContext(false);
        Configuration config = context.getConfiguration();

        // Make sure any existing loggers are removed
        for (Appender appender : context.getRootLogger().getAppenders().values()) {
            context.getRootLogger().removeAppender(appender);
        }

        // Add the appenders to the root logger
        Logger rootLogger = context.getRootLogger();
        LoggerConfig rootLoggerConfig = config.getRootLogger();
        config.addLogger(rootLogger.getName(), rootLoggerConfig);
        context.updateLoggers(config);

        PatternLayout loggingLayoutReduced = PatternLayout.newBuilder()
                .withPattern("%msg %throwable{short.message}%n")
                .build();

        // Stdout appender
        ConsoleAppender stdOut = ConsoleAppender.newBuilder()
                .setName("stdout")
                .setLayout(loggingLayoutReduced)
                .build();
        stdOut.start();

        if (logFile != null) {
            PatternLayout loggingLayoutFull = PatternLayout.newBuilder()
                    .withPattern("[%level]\t%d{ABSOLUTE} - %msg%n")
                    .build();


            // Log file appender
            FileAppender file = FileAppender.newBuilder()
                    .setName("file")
                    .setLayout(loggingLayoutFull)
                    .withFileName(logFile.getCanonicalPath())
                    .build();
            file.start();

            rootLoggerConfig.addAppender(file, loggingLevel, null);

            // Create a separator in the log, so separate runs are easier to spot
            LOGGER.info("");
            LOGGER.info("Datg converter " + VERSION);
            LOGGER.info("Current date and time: " + startDateTime);
        }

        rootLoggerConfig.addAppender(stdOut, Level.INFO, null);
    }


}
