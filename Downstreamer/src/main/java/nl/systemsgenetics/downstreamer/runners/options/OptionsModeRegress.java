package nl.systemsgenetics.downstreamer.runners.options;

import org.apache.commons.cli.*;

import java.io.File;

public class OptionsModeRegress extends OptionsBase {

    private final File responseVariable;
    private final File explanatoryVariables;
    private final File covariates;
    private final File sigma;
    private final File eigenvectors;
    private final File eigenvalues;

    private final double percentageOfVariance;

    static {
        // The response variable
        OptionBuilder.withArgName("path");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("The response variable for regression (the gene-pvalues).");
        OptionBuilder.withLongOpt("response");
        OptionBuilder.isRequired();
        OPTIONS.addOption(OptionBuilder.create("y"));

        // The explanatory variable
        OptionBuilder.withArgName("path");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("The explanatory variable for regression (the pathway). Each column in this matrix is treated as a seperate regression.");
        OptionBuilder.withLongOpt("explain");
        OptionBuilder.isRequired();
        OPTIONS.addOption(OptionBuilder.create("x"));

        // The explanatory variable
        OptionBuilder.withArgName("cov");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("The optional covariates for regression.");
        OptionBuilder.withLongOpt("cov");
        OPTIONS.addOption(OptionBuilder.create("c"));

        // The correlation matrix to use for eigen decomp
        OptionBuilder.withArgName("path");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("The correlation matrix to be used as Sigma (the gene-gene correlation matrix). If not specified, must provide -U and -L");
        OptionBuilder.withLongOpt("sigma");
        OPTIONS.addOption(OptionBuilder.create("s"));

        // The eigen decomp
        OptionBuilder.withArgName("path");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("The eigenvector matrix, assumed to be sorted on decreasing order of eigenvalues. Should have dimnames");
        OptionBuilder.withLongOpt("eigenvectors");
        OPTIONS.addOption(OptionBuilder.create("u"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("The eigenvalues, assumed to be sorted on decreasing order. Should be a 1 column matrix with dimnames.");
        OptionBuilder.withLongOpt("eigenvalues");
        OPTIONS.addOption(OptionBuilder.create("l"));

        OptionBuilder.withArgName("double between 0-1");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("The number of eigenvectors explaining the percentage of variance to use in the regression.");
        OptionBuilder.withLongOpt("perc-variance");
        OptionBuilder.isRequired();
        OPTIONS.addOption(OptionBuilder.create("p"));
    }

    public OptionsModeRegress(String... args) throws ParseException {
        super(args);

        final CommandLineParser parser = new RelaxedParser();
        final CommandLine commandLine = parser.parse(OPTIONS, args, false);

        // Check X and Y
        responseVariable = new File(commandLine.getOptionValue('y'));
        if (!responseVariable.exists()) throw new ParseException("File in -y does not exist");

        explanatoryVariables = new File(commandLine.getOptionValue('x'));
        if (!explanatoryVariables.exists()) throw new ParseException("File in -x does not exist");

        // Check for either sigma or for U+L
        if (commandLine.hasOption('s')) {
            sigma = new File(commandLine.getOptionValue('s'));
            if (!explanatoryVariables.exists()) throw new ParseException("File in -s does not exist");

            eigenvectors = null;
            eigenvalues = null;
        } else if (commandLine.hasOption('u') && commandLine.hasOption('u')) {
            eigenvectors = new File(commandLine.getOptionValue('u'));
            if (!eigenvectors.exists()) throw new ParseException("File in -u does not exist");

            eigenvalues = new File(commandLine.getOptionValue('l'));
            if (!eigenvalues.exists()) throw new ParseException("File in -l does not exist");

            sigma = null;
        } else {
            throw new ParseException("Either specify -s or -U + -L");
        }

        // Optional covariates
        if (commandLine.hasOption('c')) {
            covariates = new File(commandLine.getOptionValue('c'));
            if (!explanatoryVariables.exists()) throw new ParseException("File in -c does not exist");
        } else {
            covariates = null;
        }

        percentageOfVariance = Double.parseDouble(commandLine.getOptionValue("p"));

    }

    public File getResponseVariable() {
        return responseVariable;
    }

    public File getExplanatoryVariables() {
        return explanatoryVariables;
    }

    public File getCovariates() {
        return covariates;
    }

    public File getSigma() {
        return sigma;
    }

    public File getEigenvectors() {
        return eigenvectors;
    }

    public File getEigenvalues() {
        return eigenvalues;
    }

    public boolean hasCovariates() {
        return covariates != null;
    }

    public double getPercentageOfVariance() {
        return percentageOfVariance;
    }

    public boolean hasSigma() {
        return sigma != null;
    }
}
