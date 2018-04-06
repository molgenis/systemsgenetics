/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.deelenp.genotypeharmonizer;

import java.io.File;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.log4j.Logger;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.GenotypedDataWriterFormats;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;

/**
 *
 * @author Patrick Deelen
 */
public class GenotypeHarmonizerParamaters {

    private static final Logger LOGGER;
    private static final Options OPTIONS;
    private final String[] inputBasePaths;
    private final RandomAccessGenotypeDataReaderFormats inputType;
    private final String refBasePath;
    private final RandomAccessGenotypeDataReaderFormats refType;
    private final String outputBasePath;
    private final GenotypedDataWriterFormats outputType;
    private final boolean debugMode;
    private final boolean updateId;
    private final boolean matchRefAllele;
    private final boolean ambiguousSnpFilter;
    private final int minSnpsToAlignOn;
    private final int flankSnpsToConsider;
    private final double minLdToIncludeAlign;
    private final double maxMafForMafAlignment;
    private final double minimumPosteriorProbability;
    private final String forceSeqName;
    private final boolean ldCheck;
    private final boolean keep;
    private final File variantFilterListFile;
    private final File variantPosFilterListFile;
    private final String seqFilterIn;
    private final File sampleFilterListFile;
    private final File logFile;
    private final File snpUpdateFile;
    private final File snpLogFile;
    private final double minHwePvalue;
    private final float minCallRate;
    private final float minMAF;
    private final double minMachR2;
    /**
     * The default minimum number of SNPs that must have LD above minimum LD
     * before doing alignment based on LD
     */
    private static final int DEFAULT_MIN_VARIANTS_TO_ALIGN_ON = 3;
    /**
     * The default number of SNPs on either flank to consider for LD alignment.
     * Only SNPs with LD above minimum LD will be used
     */
    private static final int DEFAULT_FLANK_VARIANTS_TO_CONSIDER = 100;
    /**
     * The default minimum LD before using a SNP for LD alignment
     */
    private static final double DEFAULT_MIN_LD_TO_INCLUDE_ALIGN = 0.3;
    /**
     * The default maximum LD before the minor allele frequency (MAF) is used as
     * backup for alignment
     */
    private static final double DEFAULT_MAX_MAF_FOR_MAF_ALIGNMENT = 0;
    /**
     * The default minimum posterior probability to call genotypes
     */
    private static final double DEFAULT_MINIMUM_POSTERIOR_PROBABILITY = 0.4;

    static {

        LOGGER = Logger.getLogger(GenotypeHarmonizerParamaters.class);

        OPTIONS = new Options();

        Option option;

        option = OptionBuilder.withArgName("basePath")
                .hasArgs()
                .withDescription("The base path of the data to align. The extensions are determined based on the input data type.")
                .withLongOpt("input")
                .isRequired()
                .create("i");
        OPTIONS.addOption(option);

        option = OptionBuilder.withArgName("basePath")
                .hasArg()
                .withDescription("The base bath of the reference data. The extensions are determined based on the reference data type.")
                .withLongOpt("ref")
                .create("r");
        OPTIONS.addOption(option);

        option = OptionBuilder.withArgName("type")
                .hasArg()
                .withDescription("The input data type. If not defined will attempt to automatically select the first matching dataset on the specified path\n"
                        + "* PED_MAP - plink PED MAP files.\n"
                        + "* PLINK_BED - plink BED BIM FAM files.\n"
                        + "* VCF - bgziped vcf with tabix index file\n"
                        + "* VCFFOLDER - matches all bgziped vcf files + tabix index in a folder\n"
                        + "* SHAPEIT2 - shapeit2 phased haplotypes .haps & .sample\n"
                        + "* GEN - Oxford .gen & .sample\n"
                        + "* TRITYPER - TriTyper format folder")
                .withLongOpt("inputType")
                .create("I");
        OPTIONS.addOption(option);

        option = OptionBuilder.withArgName("type")
                .hasArg()
                .withDescription("The reference data type. If not defined will attempt to automatically select the first matching dataset on the specified path\n"
                        + "* PED_MAP - plink PED MAP files.\n"
                        + "* PLINK_BED - plink BED BIM FAM files.\n"
                        + "* VCF - bgziped vcf with tabix index file\n"
                        + "* VCF_FOLDER - matches all bgziped vcf files + tabix index in a folder\n"
                        + "* SHAPEIT2 - shapeit2 phased haplotypes .haps & .sample\n"
                        + "* GEN - Oxford .gen & .sample\n"
                        + "* TRITYPER - TriTyper format folder")
                .withLongOpt("refType")
                .create("R");
        OPTIONS.addOption(option);

        option = OptionBuilder.withArgName("basePath")
                .hasArg()
                .withDescription("The output bash path")
                .withLongOpt("output")
                .isRequired()
                .create("o");
        OPTIONS.addOption(option);

        option = OptionBuilder.withArgName("type")
                .hasArg()
                .withDescription("The output data type. Defaults to --inputType or to PLINK_BED if there is no writer for the impute type.\n"
                        + "* PED_MAP - plink PED MAP files\n"
                        + "* PLINK_BED - plink BED BIM FAM files\n"
                        + "* SHAPEIT2 - shapeit2 phased haplotypes\n"
                        + "* GEN - Oxford .gen & .sample\n"
                        + "* TRITYPER - TriTyper format folder")
                .withLongOpt("outputType")
                .create("O");
        OPTIONS.addOption(option);

        option = OptionBuilder.withArgName("int")
                .hasArg()
                .withDescription("Number of variants on either flank to consider using for LD-strand alignment. Must be equal or larger than --min-variants. Defaults to " + DEFAULT_FLANK_VARIANTS_TO_CONSIDER)
                .withLongOpt("variants")
                .create("v");
        OPTIONS.addOption(option);

        option = OptionBuilder.withArgName("double")
                .hasArg()
                .withDescription("Minimum LD between variant to align or check and a flanking variants in both input as reference. Defaults to " + DEFAULT_MIN_LD_TO_INCLUDE_ALIGN + ". It is NOT recommend to set this to zero")
                .withLongOpt("min-ld")
                .create("l");
        OPTIONS.addOption(option);

        option = OptionBuilder.withArgName("int")
                .hasArg()
                .withDescription("Minimum number of variants above ld-cutoff to do LD alignment. Variants that do not meet this requirement are excluded. Defaults to " + DEFAULT_MIN_VARIANTS_TO_ALIGN_ON + ". Min value: " + GenotypeHarmonizer.MIN_MIN_VARIANTS_TO_ALIGN_ON)
                .withLongOpt("min-variants")
                .create("m");
        OPTIONS.addOption(option);

        option = OptionBuilder.withArgName("boolean")
                .withDescription("Check ld structure of all variants after alignment and exclude variants that deviate. This option negates the effect of --mafAlign")
                .withLongOpt("check-ld")
                .create("c");
        OPTIONS.addOption(option);

        option = OptionBuilder.withArgName("boolean")
                .withDescription("Activate debug mode. This will result in a more verbose log file")
                .withLongOpt("debug")
                .create("d");
        OPTIONS.addOption(option);

        option = OptionBuilder.withArgName("boolean")
                .withDescription("Update variants IDs of study data to match reference data")
                .withLongOpt("update-id")
                .create("id");
        OPTIONS.addOption(option);
        
        option = OptionBuilder.withArgName("boolean")
                .withDescription("Filter out ambigous SNPs")
                .withLongOpt("ambiguousSnpFilter")
                .create("asf");
        OPTIONS.addOption(option);
        
        option = OptionBuilder.withArgName("boolean")
                .withDescription("Match reference allele to harmonizing panel")
                .withLongOpt("update-reference-allele")
                .create("ura");
        OPTIONS.addOption(option);

        option = OptionBuilder.withArgName("boolean")
                .withDescription("Keep variants not present in reference data")
                .withLongOpt("keep")
                .create("k");
        OPTIONS.addOption(option);

        option = OptionBuilder.withArgName("string")
                .hasArg()
                .withDescription("Shapeit2 does not output the sequence name in the first column of the haplotype file. Use this option to force the chromosome for all variants. This option is only valid in combination with --inputType SHAPEIT2")
                .withLongOpt("forceChr")
                .create("f");
        OPTIONS.addOption(option);

        option = OptionBuilder.withArgName("double")
                .hasArg()
                .withDescription("If there are not enough variants in LD and the minor allele frequency (MAF) of a variant <= the specified value in both study as in reference then the minor allele can be used as a backup for alignment. Defaults to " + DEFAULT_MAX_MAF_FOR_MAF_ALIGNMENT)
                .withLongOpt("mafAlign")
                .create("ma");
        OPTIONS.addOption(option);

        option = OptionBuilder.withArgName("double")
                .hasArg()
                .withDescription("The minimum posterior probability to call genotypes in the input data " + DEFAULT_MINIMUM_POSTERIOR_PROBABILITY)
                .withLongOpt("inputProb")
                .create("ip");
        OPTIONS.addOption(option);

        option = OptionBuilder.withArgName("double")
                .hasArg()
                .withDescription("The minimum minor allele frequency to include variant from input data")
                .withLongOpt("mafFilter")
                .create("mf");
        OPTIONS.addOption(option);

        option = OptionBuilder.withArgName("double")
                .hasArg()
                .withDescription("The minimum MACH R2 measure to include SNPs")
                .withLongOpt("machR2Filter")
                .create("mrf");
        OPTIONS.addOption(option);

        option = OptionBuilder.withArgName("double")
                .hasArg()
                .withDescription("The minimum call rate to include variant from input data")
                .withLongOpt("callRateFilter")
                .create("cf");
        OPTIONS.addOption(option);

        option = OptionBuilder.withArgName("double")
                .hasArg()
                .withDescription("The minimum hardy weinberg equilibrium p-value to include variant from input data")
                .withLongOpt("hweFilter")
                .create("hf");
        OPTIONS.addOption(option);

        option = OptionBuilder.withArgName("string")
                .hasArg()
                .withDescription("Filter input data on chromosome")
                .withLongOpt("chrFilter")
                .create("ch");
        OPTIONS.addOption(option);

        option = OptionBuilder.withArgName("string")
                .hasArg()
                .withDescription("Path to file with variant IDs to include from input data.")
                .withLongOpt("variantFilterList")
                .create("vf");
        OPTIONS.addOption(option);

        option = OptionBuilder.withArgName("string")
                .hasArg()
                .withDescription("Path to file with samples IDs to include from input data. For plink data only the sample id (column 2) is used")
                .withLongOpt("sampleFilterList")
                .create("sf");
        OPTIONS.addOption(option);

        option = OptionBuilder.withArgName("string")
                .hasArg()
                .withDescription("Path to file with variant CHR\tPOS or CHR:POS to include from input data.")
                .withLongOpt("variantPosFilterList")
                .create("pf");
        OPTIONS.addOption(option);

    }

    public GenotypeHarmonizerParamaters(String... args) throws ParseException {

        CommandLineParser parser = new PosixParser();
        final CommandLine commandLine = parser.parse(OPTIONS, args, false);

        inputBasePaths = commandLine.getOptionValues('i');

        try {
            if (commandLine.hasOption('I')) {
                inputType = RandomAccessGenotypeDataReaderFormats.valueOfSmart(commandLine.getOptionValue('I').toUpperCase());
            } else {
                if (inputBasePaths[0].endsWith(".vcf")) {
                    throw new ParseException("Only vcf.gz is supported. Please see manual on how to do create a vcf.gz file.");
                }
                try {
                    inputType = RandomAccessGenotypeDataReaderFormats.matchFormatToPath(inputBasePaths[0]);
                } catch (GenotypeDataException e) {
                    throw new ParseException("Unable to determine input type based on specified path. Please specify --inputType");
                }
            }
        } catch (IllegalArgumentException e) {
            throw new ParseException("Error parsing --inputType \"" + commandLine.getOptionValue('I') + "\" is not a valid input data format");
        }

        if (commandLine.hasOption('r')) {

            refBasePath = commandLine.getOptionValue('r');

            try {
                if (commandLine.hasOption('R')) {
                    refType = RandomAccessGenotypeDataReaderFormats.valueOfSmart(commandLine.getOptionValue('R').toUpperCase());
                } else {
                    if (refBasePath.endsWith(".vcf")) {
                        throw new ParseException("Only vcf.gz is supported. Please see manual on how to do create a vcf.gz file.");
                    }
                    try {
                        refType = RandomAccessGenotypeDataReaderFormats.matchFormatToPath(refBasePath);
                    } catch (GenotypeDataException e) {
                        throw new ParseException("Unable to determine reference type based on specified path. Please specify --refType");
                    }
                }

            } catch (IllegalArgumentException e) {
                throw new ParseException("Error parsing --refType \"" + commandLine.getOptionValue('R') + "\" is not a valid reference data format");
            }

        } else {
            refBasePath = null;
            refType = null;
        }

        outputBasePath = commandLine.getOptionValue('o');

        if (commandLine.hasOption('O')) {
            try {
                outputType = GenotypedDataWriterFormats.valueOf(commandLine.getOptionValue('O').toUpperCase());
            } catch (IllegalArgumentException e) {
                throw new ParseException("Error parsing --outputType \"" + commandLine.getOptionValue('O') + "\" is not a valid output data format");
            }
        } else {
            GenotypedDataWriterFormats outputTypeTmp;
            try {
                outputTypeTmp = GenotypedDataWriterFormats.valueOf(inputType.name());
            } catch (IllegalArgumentException e) {
                outputTypeTmp = GenotypedDataWriterFormats.PLINK_BED;
            }
            outputType = outputTypeTmp;
        }

        debugMode = commandLine.hasOption('d');
        updateId = commandLine.hasOption("id");
        matchRefAllele = commandLine.hasOption("ura");
        ambiguousSnpFilter = commandLine.hasOption("asf");

        try {
            minSnpsToAlignOn = commandLine.hasOption('m') ? Integer.parseInt(commandLine.getOptionValue('m')) : DEFAULT_MIN_VARIANTS_TO_ALIGN_ON;
        } catch (NumberFormatException e) {
            throw new ParseException("Error parsing --min-variants \"" + commandLine.getOptionValue('m') + "\" is not an int");
        }

        try {
            flankSnpsToConsider = commandLine.hasOption('v') ? Integer.parseInt(commandLine.getOptionValue('v')) : DEFAULT_FLANK_VARIANTS_TO_CONSIDER;
        } catch (NumberFormatException e) {
            throw new ParseException("Error parsing --variants \"" + commandLine.getOptionValue('v') + "\" is not an int");
        }

        try {
            minLdToIncludeAlign = commandLine.hasOption('l') ? Double.parseDouble(commandLine.getOptionValue('l')) : DEFAULT_MIN_LD_TO_INCLUDE_ALIGN;
        } catch (NumberFormatException e) {
            throw new ParseException("Error parsing --min-ld \"" + commandLine.getOptionValue('l') + "\" is not a double");
        }

        try {
            maxMafForMafAlignment = commandLine.hasOption("ma") ? Double.parseDouble(commandLine.getOptionValue("ma")) : DEFAULT_MAX_MAF_FOR_MAF_ALIGNMENT;
        } catch (NumberFormatException e) {
            throw new ParseException("Error parsing --mafAlign \"" + commandLine.getOptionValue("ma") + "\" is not a double");
        }

        try {
            minimumPosteriorProbability = commandLine.hasOption("ip") ? Double.parseDouble(commandLine.getOptionValue("ip")) : DEFAULT_MINIMUM_POSTERIOR_PROBABILITY;
        } catch (NumberFormatException e) {
            throw new ParseException("Error parsing --inputProb \"" + commandLine.getOptionValue("ip") + "\" is not an double");
        }

        forceSeqName = commandLine.hasOption('f') ? commandLine.getOptionValue('f') : null;
        ldCheck = commandLine.hasOption('c');
        keep = commandLine.hasOption('k');
        if (outputType == GenotypedDataWriterFormats.TRITYPER) {
            File outputFolder = new File(outputBasePath);
            logFile = new File(outputFolder, "GenotypeHarmonizer.log");
            snpUpdateFile = new File(outputFolder, "GenotypeHarmonizer_idUpdates.txt");
            snpLogFile = new File(outputFolder, "GenotypeHarmonizer_snpLog.log");
        } else {
            logFile = new File(outputBasePath + ".log");
            snpUpdateFile = new File(outputBasePath + "_idUpdates.txt");
            snpLogFile = new File(outputBasePath + "_snpLog.log");
        }

        variantFilterListFile = commandLine.hasOption("vf") ? new File(commandLine.getOptionValue("vf")) : null;
        variantPosFilterListFile = commandLine.hasOption("pf") ? new File(commandLine.getOptionValue("pf")) : null;
        seqFilterIn = commandLine.hasOption("ch") ? commandLine.getOptionValue("ch") : null;
        sampleFilterListFile = commandLine.hasOption("sf") ? new File(commandLine.getOptionValue("sf")) : null;

        try {
            minHwePvalue = commandLine.hasOption("hf") ? Double.parseDouble(commandLine.getOptionValue("hf")) : 0.0D;
        } catch (NumberFormatException e) {
            throw new ParseException(new StringBuilder().append("Error parsing --hweFilter \"").append(commandLine.getOptionValue("hw")).append("\" is not a double").toString());
        }

        try {
            minMAF = commandLine.hasOption("mf") ? Float.parseFloat(commandLine.getOptionValue("mf")) : 0.0F;
        } catch (NumberFormatException e) {
            throw new ParseException(new StringBuilder().append("Error parsing --mafFilter \"").append(commandLine.getOptionValue("mf")).append("\" is not a double").toString());
        }

        try {
            minMachR2 = commandLine.hasOption("mrf") ? Double.parseDouble(commandLine.getOptionValue("mrf")) : 0.0d;
        } catch (NumberFormatException e) {
            throw new ParseException(new StringBuilder().append("Error parsing --machR2Filter \"").append(commandLine.getOptionValue("mrf")).append("\" is not a double").toString());
        }

        try {
            minCallRate = commandLine.hasOption("cf") ? Float.parseFloat(commandLine.getOptionValue("cf")) : 0.0F;
        } catch (NumberFormatException e) {
            throw new ParseException(new StringBuilder().append("Error parsing --callRateFilter \"").append(commandLine.getOptionValue("cf")).append("\" is not an double").toString());
        }

    }

    public static void printHelp() {
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp(" ", OPTIONS);
    }

    public void printOptions() {

        StringBuilder inputPaths = new StringBuilder();
        for (String path : inputBasePaths) {
            inputPaths.append(path);
            inputPaths.append(' ');
        }

        System.out.println("Interpreted arguments: ");
        System.out.println(" - Input base path: " + inputPaths);
        LOGGER.info("Input base path: " + inputPaths);
        System.out.println(" - Input data type: " + inputType.getName());
        LOGGER.info("Input data type: " + inputType.getName());
        System.out.println(" - Output base path: " + outputBasePath);
        LOGGER.info("Output base path: " + outputBasePath);
        System.out.println(" - Output data type: " + outputType.getName());
        LOGGER.info("Output data type: " + outputType.getName());

        if (refBasePath != null) {
            System.out.println(" - Reference base path: " + refBasePath);
            LOGGER.info("Reference base path: " + refBasePath);
            System.out.println(" - Reference data type: " + refType.getName());
            LOGGER.info("Reference data type: " + refType.getName());

            System.out.println(" - Number of flank variants to consider for LD alignment: " + flankSnpsToConsider);
            LOGGER.info("Number of flank variants to consider for LD alignment: " + flankSnpsToConsider);
            System.out.println(" - Minimum LD of flanking variants before using for LD alignment: " + minLdToIncludeAlign);
            LOGGER.info("Minimum LD of flanking variants before using for LD alignment: " + minLdToIncludeAlign);
            System.out.println(" - Minimum number of variants needed to for LD alignment: " + minSnpsToAlignOn);
            LOGGER.info("Minimum number of variants needed to for LD alignment: " + minSnpsToAlignOn);
            System.out.println(" - Maximum MAF of variants to use minor allele as backup for alignment: " + maxMafForMafAlignment);
            LOGGER.info("Maximum MAF of variants to use minor allele as backup for alignment: " + maxMafForMafAlignment);
            System.out.println(" - Update study IDs: " + (updateId ? "yes" : "no"));
            LOGGER.info("Update study variant IDs: " + (updateId ? "yes" : "no"));
            System.out.println(" - Match study reference alleles: " + (matchRefAllele ? "yes" : "no"));
            LOGGER.info("Match study reference alleles: " + (matchRefAllele ? "yes" : "no"));
            System.out.println(" - Keep variants not in reference data: " + (keep ? "yes" : "no"));
            LOGGER.info("Keep variants not in reference data: " + (keep ? "yes" : "no"));
        } else {
            System.out.println(" - Reference base path not set, not performing harmonization.");
            LOGGER.info("Reference base path not set, not performing harmonization.");
        }

        System.out.println(" - Minimum posterior probability for input data: " + minimumPosteriorProbability);
        LOGGER.info("Minimum posterior probability for input data: " + minimumPosteriorProbability);

        System.out.println(" - LD checker " + (ldCheck ? "on" : "off"));
        LOGGER.info("LD checker " + (ldCheck ? "on" : "off"));

        System.out.println(" - Force input sequence name: " + (forceSeqName == null ? "not forcing" : "forcing to: " + forceSeqName));
        LOGGER.info("Force input sequence name: " + (forceSeqName == null ? "not forcing" : "forcing to: " + forceSeqName));

        if (sampleFilterListFile != null) {
            LOGGER.info("Filter input data to samples present in: " + sampleFilterListFile);
            System.out.println(" - Filter input data to samples present in: " + sampleFilterListFile);
        }

        if (variantFilterListFile != null) {
            LOGGER.info("Filter input data to variants present in: " + variantFilterListFile);
            System.out.println(" - Filter input data to variants present in: " + variantFilterListFile);
        }
        if (ambiguousSnpFilter) {
            LOGGER.info("Filter A/T & C/G SNPs.");
            System.out.println(" - Filter A/T & C/G SNPs. ");
        }
        if (variantPosFilterListFile != null) {
            LOGGER.info("Filter input data to variants with postions present in: " + variantPosFilterListFile);
            System.out.println(" - Filter input data to variants with postions present in: " + variantPosFilterListFile);
        }

        if (minMAF > 0) {
            LOGGER.info("Filter input data on minimum MAF: " + minMAF);
            System.out.println(" - Filter input data on minimum MAF: " + minMAF);
        }

        if (minHwePvalue > 0) {
            LOGGER.info("Filter input data on minimum HWE p-value: " + minHwePvalue);
            System.out.println(" - Filter input data on minimum HWE p-value: " + minHwePvalue);
        }

        if (minCallRate > 0) {
            LOGGER.info("Filter input data on minimum variant call-rate: " + minCallRate);
            System.out.println(" - Filter input data on minimum variant call-rate: " + minCallRate);
        }

        if (minMachR2 > 0) {
            LOGGER.info("Filter input data on minimum MACH R2 measure: " + minMachR2);
            System.out.println(" - Filter input data on minimum MACH R2 measure: " + minMachR2);
        }

        if (seqFilterIn != null) {
            LOGGER.info("Filter input data on variant on: " + seqFilterIn);
            System.out.println(" - Filter input data on variant on: " + seqFilterIn);
        }

        LOGGER.info("Debug mode: " + (debugMode ? "on" : "off"));

        System.out.println();

    }

    public String[] getInputBasePaths() {
        return inputBasePaths;
    }

    public String getRefBasePath() {
        return refBasePath;
    }

    public RandomAccessGenotypeDataReaderFormats getRefType() {
        return refType;
    }

    public String getOutputBasePath() {
        return outputBasePath;
    }

    public GenotypedDataWriterFormats getOutputType() {
        return outputType;
    }

    public boolean isDebugMode() {
        return debugMode;
    }

    public boolean isUpdateId() {
        return updateId;
    }

    public int getMinSnpsToAlignOn() {
        return minSnpsToAlignOn;
    }

    public int getFlankSnpsToConsider() {
        return flankSnpsToConsider;
    }

    public double getMinLdToIncludeAlign() {
        return minLdToIncludeAlign;
    }

    public double getMaxMafForMafAlignment() {
        return maxMafForMafAlignment;
    }

    public double getMinimumPosteriorProbability() {
        return minimumPosteriorProbability;
    }

    public String getForceSeqName() {
        return forceSeqName;
    }

    public boolean isLdCheck() {
        return ldCheck;
    }

    public boolean isKeep() {
        return keep;
    }

    public File getVariantFilterListFile() {
        return variantFilterListFile;
    }
    
    public boolean filterOnAmbiguousSnps(){
        return ambiguousSnpFilter;
    }

    public File getVariantPosFilterListFile() {
        return variantPosFilterListFile;
    }

    public String getSeqFilterIn() {
        return seqFilterIn;
    }

    public File getSampleFilterListFile() {
        return sampleFilterListFile;
    }

    public File getLogFile() {
        return logFile;
    }

    public File getSnpUpdateFile() {
        return snpUpdateFile;
    }

    public double getMinHwePvalue() {
        return minHwePvalue;
    }

    public float getMinCallRate() {
        return minCallRate;
    }

    public float getMinMAF() {
        return minMAF;
    }

    public RandomAccessGenotypeDataReaderFormats getInputType() {
        return inputType;
    }

    public File getSnpLogFile() {
        return snpLogFile;
    }

    public double getMinMachR2() {
        return minMachR2;
    }

    public boolean getMatchRefAllele() {
        return matchRefAllele;
    }
}
