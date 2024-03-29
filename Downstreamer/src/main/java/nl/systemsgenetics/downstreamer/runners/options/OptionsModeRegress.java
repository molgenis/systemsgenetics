package nl.systemsgenetics.downstreamer.runners.options;

import org.apache.commons.cli.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.File;

public class OptionsModeRegress extends OptionsBase {
    private static final Logger LOGGER = LogManager.getLogger(OptionsBase.class);

    private final File responseVariable;
    private final File explanatoryVariables;
    private final File covariates;
    private final File sigma;
    private final File eigenvectors;
    private final File eigenvalues;
    private final File columnIncludeFilter;
    private final File rowIncludeFilter;
    private final File genes;

    private final double percentageOfVariance;
    private final boolean useJblas;
    private final boolean fitIntercept;
    private final boolean centerAndScale;

    private final boolean inverseNormalY;
    private final boolean inverseNormalX;
    private final boolean inverseNormalC;

    private final boolean regressCovariates;

    static {
        // The response variable
        OptionBuilder.withArgName("path");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("The response variable for regression (the gene-pvalues). Currently supports only" +
                "one column. ");
        OptionBuilder.withLongOpt("response");
        OptionBuilder.isRequired();
        OPTIONS.addOption(OptionBuilder.create("y"));

        // The explanatory variable
        OptionBuilder.withArgName("path");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("The explanatory variable for regression (the pathway)." +
                " Each column in this matrix is treated as a seperate regression.");
        OptionBuilder.withLongOpt("explain");
        OptionBuilder.isRequired();
        OPTIONS.addOption(OptionBuilder.create("x"));

        // The covariates
        OptionBuilder.withArgName("cov");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("The optional covariates for regression.");
        OptionBuilder.withLongOpt("cov");
        OPTIONS.addOption(OptionBuilder.create("c"));

        // The correlation matrix to use for eigen decomp
        OptionBuilder.withArgName("path");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("The correlation matrix to be used as Sigma (the gene-gene correlation matrix). " +
                "If not specified, must provide -U and -L");
        OptionBuilder.withLongOpt("sigma");
        OPTIONS.addOption(OptionBuilder.create("s"));

        // The eigen decomp
        OptionBuilder.withArgName("path");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("The eigenvector matrix, assumed to be sorted on decreasing order of eigenvalues. " +
                "Should have dimnames");
        OptionBuilder.withLongOpt("eigenvectors");
        OPTIONS.addOption(OptionBuilder.create("u"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("The eigenvalues, assumed to be sorted on decreasing order. Should be a 1 column " +
                "matrix with dimnames.");
        OptionBuilder.withLongOpt("eigenvalues");
        OPTIONS.addOption(OptionBuilder.create("l"));

        OptionBuilder.withArgName("double between 0-1");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("The number of eigenvectors explaining the percentage of variance to use in" +
                " the regression.");
        OptionBuilder.withLongOpt("perc-variance");
        OptionBuilder.isRequired();
        OPTIONS.addOption(OptionBuilder.create("p"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Optional file with columns in X to subset before regression");
        OptionBuilder.withLongOpt("cols");
        OPTIONS.addOption(OptionBuilder.create("co"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Optional file with rows to subset before regression");
        OptionBuilder.withLongOpt("rows");
        OPTIONS.addOption(OptionBuilder.create("ro"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("File with gene info. " +
                "col1: geneName (ensg) col2: chr col3: startPos col4: stopPos col5: geneType col6: chrArm. " +
                "This is used to speed up computation by only considering correlations on the same chromosome arm. " +
                "NOTE: This assumes -u and -l have been calculated on a block diagonal. Otherwise you won't get a " +
                "speed benefit as there will be non-zero eigenvectors between off-diagonal pairs that need to be considered. " +
                "The results should still be valid, but executing time will suffer greatly.");
        OptionBuilder.withLongOpt("genes");
        OPTIONS.addOption(OptionBuilder.create("ge"));

        OptionBuilder.withArgName("boolean");
        OptionBuilder.withDescription("Use parallel colt's pure java implementation for larger matrix computations " +
                "instead of JBlas. Slower but might work if your BLAS/LAPACK is misbehaving or unavailable.");
        OptionBuilder.withLongOpt("use-colt");
        OPTIONS.addOption(OptionBuilder.create("uc"));

        OptionBuilder.withArgName("boolean");
        OptionBuilder.withDescription("Do not fit the intercept in the model.");
        OptionBuilder.withLongOpt("no-intercept");
        OPTIONS.addOption(OptionBuilder.create("ni"));

        OptionBuilder.withArgName("boolean");
        OptionBuilder.withDescription("Instead of including the covariates in the weighted model, regress them out first," +
                " run the weighted regression on these residuals. The covariates are removed in eigenvector space.");
        OptionBuilder.withLongOpt("regress-covariates");
        OPTIONS.addOption(OptionBuilder.create("rc"));

        OptionBuilder.withArgName("boolean");
        OptionBuilder.withDescription("Inverse normal transform y");
        OptionBuilder.withLongOpt("int-y");
        OPTIONS.addOption(OptionBuilder.create());

        OptionBuilder.withArgName("boolean");
        OptionBuilder.withDescription("Inverse normal transform x");
        OptionBuilder.withLongOpt("int-x");
        OPTIONS.addOption(OptionBuilder.create());

        OptionBuilder.withArgName("boolean");
        OptionBuilder.withDescription("Inverse normal transform c");
        OptionBuilder.withLongOpt("int-c");
        OPTIONS.addOption(OptionBuilder.create());

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

        percentageOfVariance = Double.parseDouble(commandLine.getOptionValue("p"));

        // Optional arguments
        if (commandLine.hasOption('c')) {
            covariates = new File(commandLine.getOptionValue('c'));
            if (!explanatoryVariables.exists()) throw new ParseException("File in -c does not exist");
        } else {
            covariates = null;
        }

        if (commandLine.hasOption("co")) {
            columnIncludeFilter = new File(commandLine.getOptionValue("co"));
            if (!columnIncludeFilter.exists()) throw new ParseException("-co does not exist");
        } else {
            columnIncludeFilter = null;
        }

        if (commandLine.hasOption("ro")) {
            rowIncludeFilter = new File(commandLine.getOptionValue("ro"));
            if (!rowIncludeFilter.exists()) throw new ParseException("-ro does not exist");
        } else {
            rowIncludeFilter = null;
        }

        if (commandLine.hasOption("ge")) {
            genes = new File(commandLine.getOptionValue("ge"));
            if (!genes.exists()) throw new ParseException("--genes does not exist");
        } else {
            genes = null;
        }

        // Boolean flags
        useJblas = !commandLine.hasOption("uc");
        fitIntercept = !commandLine.hasOption("ni");
        centerAndScale = false;

        regressCovariates = commandLine.hasOption("rc");
        if (regressCovariates && covariates == null) {
            throw new ParseException("Missing -c while -rc is specified.");
        }

        inverseNormalY = commandLine.hasOption("int-y");
        inverseNormalX = commandLine.hasOption("int-x");
        inverseNormalC = commandLine.hasOption("int-c");

        printOptions();
    }

    @Override
    public void printOptions() {
        super.printOptions();

        LOGGER.info(" * responseVariable: " + responseVariable.getPath());
        LOGGER.info(" * explanatoryVariables: " + explanatoryVariables.getPath());

        if (covariates != null) {
            LOGGER.info(" * covariates: " + covariates.getPath());
        }

        if (sigma != null) {
            LOGGER.info(" * covariates: " + sigma.getPath());
        }

        if (eigenvectors != null) {
            LOGGER.info(" * eigenvectors: " + eigenvectors.getPath());
        }

        if (eigenvalues != null) {
            LOGGER.info(" * eigenvalues: " + eigenvalues.getPath());
        }

        if (columnIncludeFilter != null) {
            LOGGER.info(" * columnIncludeFilter: " + columnIncludeFilter.getPath());
        }

        if (rowIncludeFilter != null) {
            LOGGER.info(" * rowIncludeFilter: " + rowIncludeFilter.getPath());
        }

        if (genes != null) {
            LOGGER.info(" * genes: " + genes.getPath());
        }

        LOGGER.info(" * percentageOfVariance: " + percentageOfVariance);
        LOGGER.info(" * useJblas: " + useJblas);
        LOGGER.info(" * fitIntercept: " + fitIntercept);
        LOGGER.info(" * regressCovariates: " + regressCovariates);

        LOGGER.info(" * centerAndScale: " + centerAndScale);
        LOGGER.info(" * inverseNormalY: " + inverseNormalY);
        LOGGER.info(" * inverseNormalX: " + inverseNormalX);
        LOGGER.info(" * inverseNormalC: " + inverseNormalC);

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

    public File getColumnIncludeFilter() {
        return columnIncludeFilter;
    }

    public File getRowIncludeFilter() {
        return rowIncludeFilter;
    }

    public File getGenes() {
        return genes;
    }

    public boolean useJblas() {
        return useJblas;
    }

    public boolean fitIntercept() {
        return fitIntercept;
    }

    public boolean centerAndScale() {
        return centerAndScale;
    }

    public boolean isInverseNormalY() {
        return inverseNormalY;
    }

    public boolean isInverseNormalX() {
        return inverseNormalX;
    }

    public boolean isInverseNormalC() {
        return inverseNormalC;
    }

    public boolean regressCovariates() {
        return regressCovariates;
    }

    public boolean hasSigma() {
        return sigma != null;
    }

    public boolean hasColumnIncludeFilter() {
        return columnIncludeFilter != null;
    }

    public boolean hasRowIncludeFilter() {
        return rowIncludeFilter != null;
    }

    public boolean hasGenes() {
        return genes != null;
    }

}
