package nl.systemsgenetics.downstreamer.runners.options;

import org.apache.commons.cli.*;

import java.io.File;

import static nl.systemsgenetics.downstreamer.runners.options.DownstreamerMode.*;

public class OptionsModeCoreg extends OptionsBase {

    // TODO: refactor to input file
    private final File gwasZscoreMatrixPath;
    private final File geneInfoFile;

    private final String[] columnsToExtract;

    private final int numberSamplesUsedForCor;
    private final int cisWindow;

    private final boolean trimGeneNames;
    private final boolean normalizeEigenvectors;
    private final boolean convertRToZscore;

    static {
        // The response variable
        OptionBuilder.withArgName("path");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("The input file, or prefix to input file, to convert. Does not perse have to be a GWAS");
        OptionBuilder.withLongOpt("gwas");
        OptionBuilder.isRequired();
        OPTIONS.addOption(OptionBuilder.create("g"));

        OptionBuilder.withArgName("boolean");
        OptionBuilder.withDescription("For mode X | X | X  trim train .## from ENSEMBL gene ids");
        OptionBuilder.withLongOpt("trimGeneNames");
        OPTIONS.addOption(OptionBuilder.create("tgn"));

        OptionBuilder.withArgName("strings");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("Column names (seperated by space) to extract when running --mode X or X");
        OptionBuilder.withLongOpt("columnsToExtract");
        OPTIONS.addOption(OptionBuilder.create("cte"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("File with gene info. col1: geneName (ensg), col2: chr, col3: startPos, col4: stopPos, col5: geneType, col6: karyotype");
        OptionBuilder.withLongOpt("genes");
        OPTIONS.addOption(OptionBuilder.create("ge"));

        // TODO: check if row normalization is valid. I agree with the per col, but why per row I don't get.
        OptionBuilder.withArgName("boolean");
        OptionBuilder.withDescription("For mode COREG_CORRELATE_GENES: First normalize the eigen vectors by row, then by col");
        OptionBuilder.withLongOpt("normalizeEigenvectors");
        OPTIONS.addOption(OptionBuilder.create("ne"));

        OptionBuilder.withArgName("boolean");
        OptionBuilder.withDescription("For mode COREG_CORRELATE_GENES: Save results as Z-scores instead of r's.");
        OptionBuilder.withLongOpt("corZscore");
        OPTIONS.addOption(OptionBuilder.create("cz"));

        OptionBuilder.withArgName("int");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("For mode R_TO_ZSCORE | ZSCORE_TO_R | INVESTIGATE_NETWORK: Specify the number of samples used to create the correlation matrix");
        OptionBuilder.withLongOpt("numberSamplesUsedForCor");
        OPTIONS.addOption(OptionBuilder.create("ns"));

        OptionBuilder.withArgName("int");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("For mode COREG_REMOVE_CIS_COEXP: Cis window is defined as [start - cwe, end + cwe]. Defaults to: 250000");
        OptionBuilder.withLongOpt("cisWindowExtent");
        OPTIONS.addOption(OptionBuilder.create("cwe"));

    }

    public OptionsModeCoreg(String[] args) throws ParseException {
        super(args);

        // Parse arguments
        final CommandLineParser parser = new PosixParser();
        final CommandLine commandLine = parser.parse(OPTIONS, args, false);

        // Required & boolean arguments
        gwasZscoreMatrixPath = new File(commandLine.getOptionValue('g'));
        trimGeneNames = commandLine.hasOption("tgn");
        normalizeEigenvectors = commandLine.hasOption("ne");
        convertRToZscore = commandLine.hasOption("cz");

        if (commandLine.hasOption("ge")) {
            geneInfoFile = new File(commandLine.getOptionValue("ge"));
        } else {
            geneInfoFile = null;
        }

        if (commandLine.hasOption("cte")) {
            columnsToExtract = commandLine.getOptionValues("cte");
        } else {
            columnsToExtract = null;
        }


        try {

            if (getMode() == COREG_RTOZSCORE || getMode() == COREG_INVESTIGATE_NETWORK) {
                if (!commandLine.hasOption("ns")) {
                    throw new ParseException("--numberSamplesUsedForCor not specified");
                } else {
                    numberSamplesUsedForCor = Integer.parseInt(commandLine.getOptionValue("ns"));
                }
            } else {
                numberSamplesUsedForCor = 0;
            }

            if (getMode() == COREG_REMOVE_CIS_COEXP) {
                cisWindow = Integer.parseInt(commandLine.getOptionValue("cwe", "250000"));
            } else {
                cisWindow = 0;
            }

        } catch (NumberFormatException e) {
            throw new ParseException("Error parsing -ns | -cwe not an int");
        }

    }

    public String getGwasZscoreMatrixPath() {
        return gwasZscoreMatrixPath.getPath();
    }

    public String[] getColumnsToExtract() {
        return columnsToExtract;
    }

    public File getGeneInfoFile() {
        return geneInfoFile;
    }

    public int getNumberSamplesUsedForCor() {
        return numberSamplesUsedForCor;
    }

    public int getCisWindow() {
        return cisWindow;
    }

    public boolean isTrimGeneNames() {
        return trimGeneNames;
    }

    public boolean isNormalizeEigenvectors() {
        return normalizeEigenvectors;
    }

    public boolean isConvertRToZscore() {
        return convertRToZscore;
    }


}
