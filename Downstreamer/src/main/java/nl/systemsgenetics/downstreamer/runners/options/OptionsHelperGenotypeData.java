package nl.systemsgenetics.downstreamer.runners.options;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;

import java.io.File;

/**
 * Helper that handles parsing of genotype data related options. Can be easily used to add these options if need be.
 * See @OptionsModeConvert for an example of this.
 */
public class OptionsHelperGenotypeData implements GenotypeFileProvider {

    protected static final Options OPTIONS;

    // Genotype options
    private final String[] genotypeBasePath;
    private final RandomAccessGenotypeDataReaderFormats genotypeType;
    private final File genotypeSamplesFile;
    private final double mafFilter;

    static {
        OPTIONS = new Options();

        // The basename of genotype data
        OptionBuilder.withArgName("basePath");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("The reference genotypes");
        OptionBuilder.withLongOpt("refGenotype");
        OPTIONS.addOption(OptionBuilder.create("r"));

        // The type of genotype data
        OptionBuilder.withArgName("basePath");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("The type of reference genotype data. Automatically chosen by default.");
        OptionBuilder.withLongOpt("refGenotypeType");
        OPTIONS.addOption(OptionBuilder.create("R"));

        // Sample filter
        OptionBuilder.withArgName("path");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Samples to include from reference genotypes");
        OptionBuilder.withLongOpt("refSamples");
        OPTIONS.addOption(OptionBuilder.create("rs"));

        // MAF filter
        OptionBuilder.withArgName("double");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Minimum MAF");
        OptionBuilder.withLongOpt("maf");
        OPTIONS.addOption(OptionBuilder.create("maf"));
    }


    public OptionsHelperGenotypeData(CommandLine commandLine) throws ParseException {

        // To avoid nullpointer and help debugging
        if (!commandLine.hasOption("r")) {
            throw new ParseException("-r must be specific when providing genotype data");
        }

        genotypeBasePath = commandLine.getOptionValues('r');

        // Sample filter
        if (commandLine.hasOption("rs")) {
            genotypeSamplesFile = new File(commandLine.getOptionValue("rs"));
        } else {
            genotypeSamplesFile = null;
        }

        // Genotype data
        try {

            if (commandLine.hasOption('R')) {
                genotypeType = RandomAccessGenotypeDataReaderFormats.valueOfSmart(commandLine.getOptionValue('R').toUpperCase());
            } else {
                if (genotypeBasePath[0].endsWith(".vcf")) {
                    throw new ParseException("Only vcf.gz is supported. Please see manual on how to do create a vcf.gz file.");
                }
                try {
                    genotypeType = RandomAccessGenotypeDataReaderFormats.matchFormatToPath(genotypeBasePath);
                } catch (GenotypeDataException e) {
                    throw new ParseException("Unable to determine reference type based on specified path. Please specify --referenceGenotypesType");
                }
            }

        } catch (IllegalArgumentException e) {
            throw new ParseException("Error parsing --referenceGenotypesType \"" + commandLine.getOptionValue('R') + "\" is not a valid reference data format");
        }

        // MAF filter
        if (commandLine.hasOption("maf")) {
            try {
                mafFilter = Double.parseDouble(commandLine.getOptionValue("maf"));
            } catch (NumberFormatException e) {
                throw new ParseException("Error parsing --maf \"" + commandLine.getOptionValue("maf") + "\" is not an double");
            }
        } else {
            mafFilter = 0;
        }
    }

    public static Options getOptions() {
        return OPTIONS;
    }

    @Override
    public String[] getGenotypeBasePath() {
        return genotypeBasePath;
    }

    @Override
    public RandomAccessGenotypeDataReaderFormats getGenotypeType() {
        return genotypeType;
    }

    @Override
    public File getGenotypeSamplesFile() {
        return genotypeSamplesFile;
    }

    @Override
    public double getMafFilter() {
        return mafFilter;
    }



}
