package nl.systemsgenetics.downstreamer.runners.options;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.ParseException;

import java.io.File;

public class OptionsTesting extends OptionsBase{

    private final File sigma;
    private final File index;
    private final boolean useJblas;

    static {

        // The correlation matrix to use for eigen decomp
        OptionBuilder.withArgName("path");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("The correlation matrix to be used as Sigma (the gene-gene correlation matrix)");
        OptionBuilder.withLongOpt("sigma");
        OPTIONS.addOption(OptionBuilder.create("s"));


        OptionBuilder.withArgName("path");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("file indicating the block diagonality");
        OptionBuilder.withLongOpt("blockindices");
        OPTIONS.addOption(OptionBuilder.create("b"));

        OptionBuilder.withArgName("boolean");
        OptionBuilder.withDescription("Use parallel colt's pure java implementation for larger matrix computations " +
                "instead of JBlas. Slower but might work if your BLAS/LAPACK is misbehaving or unavailable.");
        OptionBuilder.withLongOpt("use-colt");
        OPTIONS.addOption(OptionBuilder.create("uc"));

    }

    public OptionsTesting(String[] args) throws ParseException {
        super(args);

        final CommandLineParser parser = new RelaxedParser();
        final CommandLine commandLine = parser.parse(OPTIONS, args, false);

        sigma = new File(commandLine.getOptionValue('s'));
        index = new File(commandLine.getOptionValue('b'));
        useJblas = !commandLine.hasOption("uc");

    }

    public File getIndex() {
        return index;
    }

    public File getSigma() {
        return sigma;
    }

    public boolean useJblas() {
        return useJblas;
    }
}
