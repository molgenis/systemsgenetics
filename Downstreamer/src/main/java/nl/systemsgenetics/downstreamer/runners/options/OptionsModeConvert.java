package nl.systemsgenetics.downstreamer.runners.options;

import org.apache.commons.cli.*;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;

import java.io.File;
import java.util.Collection;

public class OptionsModeConvert extends OptionsBase implements GenotypeFileProvider {

    // TODO: refactor to input file
    private final File gwasZscoreMatrixPath;
    private final File conversionColumnIncludeFilter;
    private final File conversionRowIncludeFilter;

    private final String[] columnsToExtract;

    private final boolean trimGeneNames;
    private final boolean pvalueToZscore;

    private final OptionsHelperGenotypeData optionsHelperGenotypeData;

    static {
        // The response variable
        OptionBuilder.withArgName("path");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("The input file, or prefix to input file, to convert. Does not perse have to be a GWAS");
        OptionBuilder.withLongOpt("gwas");
        OptionBuilder.isRequired();
        OPTIONS.addOption(OptionBuilder.create("g"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Optional file with columns to select during conversion");
        OptionBuilder.withLongOpt("cols");
        OPTIONS.addOption(OptionBuilder.create("co"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Optional file with rows to select during conversion");
        OptionBuilder.withLongOpt("rows");
        OPTIONS.addOption(OptionBuilder.create("ro"));

        OptionBuilder.withArgName("boolean");
        OptionBuilder.withDescription("For mode CONVERT_EXP | CONVERT_TXT  trim train .## from ENSEMBL gene ids");
        OptionBuilder.withLongOpt("trimGeneNames");
        OPTIONS.addOption(OptionBuilder.create("tgn"));

        OptionBuilder.withArgName("boolean");
        OptionBuilder.withDescription("For mode PREPARE_GWAS | PREPARE_GWAS_MERGE convert p-values to z-scores");
        OptionBuilder.withLongOpt("pvalueToZscore");
        OPTIONS.addOption(OptionBuilder.create("p2z"));

        OptionBuilder.withArgName("strings");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("Column names (seperated by space) to extract when running --mode CONVERT_BIN or CONVERT_EXP");
        OptionBuilder.withLongOpt("columnsToExtract");
        OPTIONS.addOption(OptionBuilder.create("cte"));

        // Add all options from GenotypeDataFile
        for(Option cur : (Collection<Option>) OptionsHelperGenotypeData.getOptions().getOptions()) {
            OPTIONS.addOption(cur);
        }

    }

    public OptionsModeConvert(String[] args) throws ParseException {
        super(args);

        // Parse arguments
        final CommandLineParser parser = new PosixParser();
        final CommandLine commandLine = parser.parse(OPTIONS, args, false);

        // Required arguments
        gwasZscoreMatrixPath = new File(commandLine.getOptionValue('g'));

        // Optional arguments
        if (commandLine.hasOption("co")) {
            conversionColumnIncludeFilter = new File(commandLine.getOptionValue("co"));
        } else {
            conversionColumnIncludeFilter = null;
        }

        if (commandLine.hasOption("ro")) {
            conversionRowIncludeFilter = new File(commandLine.getOptionValue("ro"));
        } else {
            conversionRowIncludeFilter = null;
        }

        trimGeneNames = commandLine.hasOption("tgn");

        if (getMode() == DownstreamerMode.CONVERT_TXT && commandLine.hasOption("p2z")) {
            throw new RuntimeException("Use mode PREPARE_GWAS to convert p-values to z-score. (this changed in version 1.32)");
        } else {
            pvalueToZscore = commandLine.hasOption("p2z");
        }

        if (commandLine.hasOption("cte")) {
            columnsToExtract = commandLine.getOptionValues("cte");
        } else {
            columnsToExtract = null;
        }

        if (commandLine.hasOption("r")) {
            optionsHelperGenotypeData = new OptionsHelperGenotypeData(commandLine);
        } else {
            optionsHelperGenotypeData = null;
        }

    }

    public String getGwasZscoreMatrixPath() {
        return gwasZscoreMatrixPath.getPath();
    }

    public File getConversionColumnIncludeFilter() {
        return conversionColumnIncludeFilter;
    }

    public File getConversionRowIncludeFilter() {
        return conversionRowIncludeFilter;
    }

    public boolean isTrimGeneNames() {
        return trimGeneNames;
    }

    public boolean isPvalueToZscore() {
        return pvalueToZscore;
    }

    public String[] getColumnsToExtract() {
        return columnsToExtract;
    }

    @Override
    public String[] getGenotypeBasePath() {
        return optionsHelperGenotypeData.getGenotypeBasePath();
    }

    @Override
    public RandomAccessGenotypeDataReaderFormats getGenotypeType() {
        return optionsHelperGenotypeData.getGenotypeType();
    }

    @Override
    public File getGenotypeSamplesFile() {
        return optionsHelperGenotypeData.getGenotypeSamplesFile();
    }

    @Override
    public double getMafFilter() {
        return optionsHelperGenotypeData.getMafFilter();
    }
}
