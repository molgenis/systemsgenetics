/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.bondermj.pcoa;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.log4j.Logger;

/**
 *
 * @author Patrick Deelen
 */
public class PcoaParamaters {

    public enum MatrixType {
        CORRELATION, COVARIATION, BRAYCURTIS, CITYBLOCK, NULL
    }

    public enum PackageToUse {
        COLT, MTJ, NULL
    }

    private static final Logger LOGGER;
    private static final Options OPTIONS;
    private final String fileIn;
    private final String numberComponentsToCalc;
    private final boolean rank;
    private final boolean overRows;
    private final boolean center;
    private final boolean scale;
    private final MatrixType matrixType;
    private final PackageToUse packageToUse;
    private boolean multiThreading = false;
    private final int threads;

    private static final String DEFAULT_numberComponents = "all";
    private static final MatrixType DEFAULT_analysisType = MatrixType.CORRELATION;
    private static final PackageToUse DEFAULT_package = PackageToUse.COLT;
    private static final boolean DEFAULT_overRows = false;
    private static final int DEFAULT_nrThreads = 1;

    static {

        LOGGER = Logger.getLogger(PcoaParamaters.class);

        OPTIONS = new Options();

        Option option;

        option = OptionBuilder.withArgName("inputFile")
                .hasArgs()
                .withDescription("The file to conduct the PCA on.")
                .withLongOpt("input")
                .isRequired()
                .create("i");
        OPTIONS.addOption(option);

        option = OptionBuilder.withArgName("String")
                .hasArg()
                .withDescription("Number components to calculate. Defaults to " + DEFAULT_numberComponents)
                .withLongOpt("nrComponents")
                .create("nc");
        OPTIONS.addOption(option);

        option = OptionBuilder.withArgName("boolean")
                .withDescription("PCA over rows. Defaults to " + DEFAULT_overRows + ". If false it will be conducted over colums (default).")
                .withLongOpt("overRows")
                .create("r");
        OPTIONS.addOption(option);
        
        option = OptionBuilder.withArgName("boolean")
                .withDescription("Rank data.")
                .withLongOpt("rank")
                .create("ra");
        OPTIONS.addOption(option);

        option = OptionBuilder.withArgName("boolean")
                .withDescription("Center input data. Defaults to " + false + ".")
                .withLongOpt("center")
                .create("c");
        OPTIONS.addOption(option);

        option = OptionBuilder.withArgName("boolean")
                .withDescription("Scale input data. Defaults to " + false + ".")
                .withLongOpt("scale")
                .create("s");
        OPTIONS.addOption(option);

        option = OptionBuilder.withArgName("String")
                .hasArg()
                .withDescription("Matrix type to perform PCA over, correlation, covariation, bray curtis or city block distances. Defaults to correlation")
                .withLongOpt("matrixType")
                .create("m");
        OPTIONS.addOption(option);

        option = OptionBuilder.withArgName("String")
                .hasArg()
                .withDescription("Package to use the eigen value decomposition from")
                .withLongOpt("package")
                .create("p");
        OPTIONS.addOption(option);

        option = OptionBuilder.withArgName("int")
                .hasArg()
                .withDescription("Number of threads to use. Defaults to: " + DEFAULT_nrThreads)
                .withLongOpt("threads")
                .create("t");
        OPTIONS.addOption(option);

    }

    public PcoaParamaters(String... args) throws ParseException {

        CommandLineParser parser = new PosixParser();
        final CommandLine commandLine = parser.parse(OPTIONS, args, false);

        fileIn = commandLine.hasOption('i') ? commandLine.getOptionValue('i') : null;
        numberComponentsToCalc = commandLine.hasOption("nc") ? commandLine.getOptionValue("nc") : DEFAULT_numberComponents;

        if (commandLine.hasOption('m')) {
            String tmp = commandLine.getOptionValue('m');
            tmp = tmp.toLowerCase();
            if (tmp.equals("correlation")) {
                matrixType = MatrixType.CORRELATION;
            } else if (tmp.equals("covariation")) {
                matrixType = MatrixType.COVARIATION;
            } else if (tmp.equals("braycurtis") || tmp.equals("bray-curtis") || tmp.equals("bray curtis")) {
                matrixType = MatrixType.BRAYCURTIS;
            } else if (tmp.equals("cityblock") || tmp.equals("city-block") || tmp.equals("city block") || tmp.equals("mannhattan distance") || tmp.equals("mannhattan-distance") || tmp.equals("mannhattandistance") || tmp.equals("L1")) {
                matrixType = MatrixType.CITYBLOCK;
            } else {
                matrixType = MatrixType.NULL;
                throw new ParseException("In valid option for -m is given: " + tmp);
            }
        } else {
            matrixType = DEFAULT_analysisType;
        }

        if (commandLine.hasOption('p')) {
            String tmp = commandLine.getOptionValue('p');
            tmp = tmp.toLowerCase();
            if (tmp.equals("colt") || tmp.equals("parallel-colt") || tmp.equals("parallel colt")) {
                packageToUse = PackageToUse.COLT;
            } else if (tmp.equals("mtj")) {
                packageToUse = PackageToUse.MTJ;
            } else {
                packageToUse = PackageToUse.NULL;
                throw new ParseException("In valid option for -p is given: " + tmp);
            }
        } else {
            packageToUse = DEFAULT_package;
        }

        if (commandLine.hasOption('c')) {
            center = true;
        } else {
            center = false;
        }
        if (commandLine.hasOption('r')) {
            overRows = true;
        } else {
            overRows = DEFAULT_overRows;
        }
        if (commandLine.hasOption("ra")) {
            rank = true;
        } else {
            rank = false;
        }
        
        if (commandLine.hasOption('s')) {
            scale = true;
        } else {
            scale = false;
        }

        try {
            threads = commandLine.hasOption('t') ? Integer.parseInt(commandLine.getOptionValue('t')) : DEFAULT_nrThreads;
            if (threads > 1) {
                multiThreading = true;
            } else {
                multiThreading = false;
            }
        } catch (NumberFormatException e) {
            throw new ParseException("Error parsing --threads \"" + commandLine.getOptionValue('t') + "\" is not an int");
        }

        try {

        } catch (NumberFormatException e) {
            throw new ParseException("Error parsing --variants \"" + commandLine.getOptionValue('v') + "\" is not an int");
        }

    }

    public static void printHelp() {
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp(" ", OPTIONS);
    }

    public static Logger getLOGGER() {
        return LOGGER;
    }

    public static Options getOPTIONS() {
        return OPTIONS;
    }

    public String getNumberComponentsToCalc() {
        return numberComponentsToCalc;
    }

    public boolean isOverRows() {
        return overRows;
    }
    
    public boolean isRank() {
        return rank;
    }
    
    public boolean isCenter() {
        return center;
    }

    public boolean isScale() {
        return scale;
    }

    public boolean isMultiThreading() {
        return multiThreading;
    }

    public int getThreads() {
        return threads;
    }

    public String getFileIn() {
        return fileIn;
    }

    public MatrixType getMatrixType() {
        return matrixType;
    }

    public PackageToUse getPackageToUse() {
        return packageToUse;
    }

}
