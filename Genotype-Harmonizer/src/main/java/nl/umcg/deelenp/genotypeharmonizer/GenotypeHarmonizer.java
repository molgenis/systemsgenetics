package nl.umcg.deelenp.genotypeharmonizer;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.DateFormat;
import java.text.NumberFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashSet;
import java.util.ResourceBundle;
import java.util.regex.Pattern;
import org.apache.commons.cli.ParseException;
import org.apache.log4j.FileAppender;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.SimpleLayout;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.GenotypeWriter;
import org.molgenis.genotype.GenotypedDataWriterFormats;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.modifiable.ModifiableGenotypeData;
import org.molgenis.genotype.multipart.IncompatibleMultiPartGenotypeDataException;
import org.molgenis.genotype.sampleFilter.SampleIdIncludeFilter;
import org.molgenis.genotype.tabix.TabixFileNotFoundException;
import org.molgenis.genotype.util.LdCalculatorException;
import org.molgenis.genotype.variantFilter.VariantCombinedFilter;
import org.molgenis.genotype.variantFilter.VariantFilterAmbigousSnp;
import org.molgenis.genotype.variantFilter.VariantFilterMachR2;
import org.molgenis.genotype.variantFilter.VariantFilterSeq;
import org.molgenis.genotype.variantFilter.VariantFilterSeqPos;
import org.molgenis.genotype.variantFilter.VariantIdIncludeFilter;
import org.molgenis.genotype.variantFilter.VariantQcChecker;

/**
 *
 * @author Patrick Deelen
 */
@SuppressWarnings("static-access")
class GenotypeHarmonizer {

    private static final String VERSION = ResourceBundle.getBundle("verion").getString("application.version");
    private static final String HEADER
            = "  /---------------------------------------\\\n"
            + "  |          Genotype Harmonizer          |\n"
            + "  |                                       |\n"
            + "  |             Patrick Deelen            |\n"
            + "  |        patrickdeelen@gmail.com        |\n"
            + "  |                                       |\n"
            + "  | Harm-Jan Westra, Joeri van der Velde, |\n"
            + "  |    Marc Jan Bonder, Erwin Winder,     |\n"
            + "  |           Dennis Hendriksen           |\n"
            + "  |      Lude Franke, Morris Swertz       |\n"
            + "  |                                       |\n"
            + "  |     Genomics Coordination Center      |\n"
            + "  |  University Medical Center Groningen  |\n"
            + "  \\---------------------------------------/";
    private static final Logger LOGGER = Logger.getLogger(GenotypeHarmonizer.class);
    private static final DateFormat DATE_TIME_FORMAT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
    /**
     * The lowest allowed minimum for the number of SNPs needed to align on
     */
    protected static final int MIN_MIN_VARIANTS_TO_ALIGN_ON = 1;
    private static final Pattern CHR_POS_SPLITTER = Pattern.compile("\\s+|:");
    public static final NumberFormat DEFAULT_NUMBER_FORMATTER = NumberFormat.getInstance();

    /**
     * @param args
     * @throws InterruptedException
     * @throws UserFriendlyException
     */
    public static void main(String... args) throws InterruptedException {

        System.out.println(HEADER);
        System.out.println();
        System.out.println("          --- Version: " + VERSION + " ---");
        System.out.println();
        System.out.println("More information: http://molgenis.org/systemsgenetics");
        System.out.println();

        Date currentDataTime = new Date();
        System.out.println("Current date and time: " + DATE_TIME_FORMAT.format(currentDataTime));
        System.out.println();

        System.out.flush(); //flush to make sure header is before errors
        Thread.sleep(25); //Allows flush to complete

        if (args.length == 0) {
            GenotypeHarmonizerParamaters.printHelp();
            System.exit(1);
        }

        final GenotypeHarmonizerParamaters paramaters;
        try {
            paramaters = new GenotypeHarmonizerParamaters(args);
        } catch (ParseException ex) {
            System.err.println("Invalid command line arguments: ");
            System.err.println(ex.getMessage());
            System.err.println();
            GenotypeHarmonizerParamaters.printHelp();
            System.exit(1);
            return;
        }

        if (paramaters.getLogFile().getParentFile() != null && !paramaters.getLogFile().getParentFile().isDirectory()) {
            if (!paramaters.getLogFile().getParentFile().mkdirs()) {
                System.err.println("Failed to create output folder: " + paramaters.getLogFile().getParent());
                System.exit(1);
            }
        }

        try {
            FileAppender logAppender = new FileAppender(new SimpleLayout(), paramaters.getLogFile().getCanonicalPath(), false);
            LOGGER.getRootLogger().removeAllAppenders();
            LOGGER.getRootLogger().addAppender(logAppender);
            if (paramaters.isDebugMode()) {
                LOGGER.setLevel(Level.DEBUG);
            } else {
                LOGGER.setLevel(Level.INFO);
            }
        } catch (IOException e) {
            System.err.println("Failed to create logger: " + e.getMessage());
            System.exit(1);
        }

        LOGGER.info(
                "\n" + HEADER);
        LOGGER.info("Version: " + VERSION);
        LOGGER.info("Current date and time: " + DATE_TIME_FORMAT.format(currentDataTime));
        LOGGER.info("Log level: " + LOGGER.getLevel());

        System.out.println(
                "Started logging");
        System.out.println();

        paramaters.printOptions();

        if (paramaters.getMinSnpsToAlignOn() < MIN_MIN_VARIANTS_TO_ALIGN_ON) {
            LOGGER.fatal("the specified --min-variants < " + MIN_MIN_VARIANTS_TO_ALIGN_ON);
            System.err.println("the specified --min-variants < " + MIN_MIN_VARIANTS_TO_ALIGN_ON);
            System.exit(1);
            return;
        }
        if (paramaters.getFlankSnpsToConsider() < paramaters.getMinLdToIncludeAlign()) {
            LOGGER.fatal("--variants < --min-variants");
            System.err.println("--variants < --min-variants");
            System.exit(1);
            return;
        }

        if (paramaters.getInputBasePaths()[0].equals(paramaters.getRefBasePath())) {
            LOGGER.fatal("Study data and reference data cannot be the same data");
            System.err.println("Study data and reference data cannot be the same data");
            System.exit(1);
            return;
        }

        if (paramaters.getInputBasePaths()[0].equals(paramaters.getOutputBasePath())) {
            LOGGER.fatal("Study input can not be the same as output");
            System.err.println("Study input can not be the same as output");
            System.exit(1);
            return;
        }

        if (paramaters.getMinCallRate() > 1) {
            LOGGER.fatal("The specified call rate is > 1. This will result in exclusion of all variants.");
            System.err.println("The specified call rate is > 1. This will result in exclusion of all variants.");
            System.exit(1);
            return;
        }

        if (paramaters.getForceSeqName() != null && paramaters.getInputType() != RandomAccessGenotypeDataReaderFormats.SHAPEIT2 && paramaters.getInputType() != RandomAccessGenotypeDataReaderFormats.GEN) {
            System.err.println("Error cannot force sequence name of: " + paramaters.getInputType().getName());
            System.exit(1);
            return;
        }
        int genotypeDataCache = paramaters.getFlankSnpsToConsider() * 4;

        VariantCombinedFilter varFilter = null;

        if ((paramaters.getMinCallRate() != 0.0F) || (paramaters.getMinMAF() != 0.0F) || (paramaters.getMinHwePvalue() != 0.0D) || (paramaters.getVariantFilterListFile() != null) || (paramaters.getSeqFilterIn() != null) || (paramaters.getVariantPosFilterListFile() != null) || paramaters.getMinMachR2() != 0.0d || paramaters.filterOnAmbiguousSnps() ) {
            varFilter = new VariantCombinedFilter();
            if (paramaters.getVariantFilterListFile() != null) {
                try {
                    HashSet<String> snps = new HashSet<String>();

                    BufferedReader variantIdFilterReader = new BufferedReader(new FileReader(paramaters.getVariantFilterListFile()));
                    String line;
                    while ((line = variantIdFilterReader.readLine()) != null) {
                        snps.add(line);
                    }
                    VariantIdIncludeFilter snpIdFilter = new VariantIdIncludeFilter(snps);
                    varFilter.add(snpIdFilter);
                } catch (FileNotFoundException ex) {
                    System.err.println("Unable to find file with variants to filter on at: " + paramaters.getVariantFilterListFile().getAbsolutePath());
                    LOGGER.fatal("Unable to find file with variants to filter on at: " + paramaters.getVariantFilterListFile().getAbsolutePath());
                    System.exit(1);
                } catch (IOException e) {
                    System.err.println("Error reading file with variants to filter on at: " + paramaters.getVariantFilterListFile().getAbsolutePath() + " error: " + e.getMessage());
                    LOGGER.fatal("Error reading file with variants to filter on at: " + paramaters.getVariantFilterListFile(), e);
                    System.exit(1);
                }
            }
            if (paramaters.getVariantPosFilterListFile() != null) {
                VariantFilterSeqPos varPosFilter = new VariantFilterSeqPos();
                try {
                    BufferedReader variantIdFilterReader = new BufferedReader(new FileReader(paramaters.getVariantPosFilterListFile()));
                    String line;
                    while ((line = variantIdFilterReader.readLine()) != null) {

                        String[] elements = CHR_POS_SPLITTER.split(line);

                        if (elements.length != 2) {
                            System.err.println("Error parsing chr pos for line: " + line + " skipping line");
                            LOGGER.error("Error parsing chr pos for line: " + line + " skipping line");
                            continue;
                        }

                        varPosFilter.addSeqPos(elements[0], Integer.parseInt(elements[1]));
                    }
                    varFilter.add(varPosFilter);

                } catch (FileNotFoundException ex) {
                    System.err.println("Unable to find file with variant positions to filter on at: " + paramaters.getVariantFilterListFile().getAbsolutePath());
                    LOGGER.fatal("Unable to find file with variant positions to filter on at: " + paramaters.getVariantPosFilterListFile());
                    System.exit(1);
                } catch (IOException e) {
                    System.err.println("Error reading file with variant positions to filter on at: " + paramaters.getVariantPosFilterListFile() + " error: " + e.getMessage());
                    LOGGER.fatal("Error reading file with variant positions to filter on at: " + paramaters.getVariantPosFilterListFile(), e);
                    System.exit(1);
                }
            }
            if ((paramaters.getMinCallRate() != 0.0F) || (paramaters.getMinMAF() != 0.0F) || (paramaters.getMinHwePvalue() != 0.0D)) {
                VariantQcChecker snpQcFilter = new VariantQcChecker(paramaters.getMinMAF(), paramaters.getMinCallRate(), paramaters.getMinHwePvalue());
                varFilter.add(snpQcFilter);
            }
            if (paramaters.getSeqFilterIn() != null) {
                VariantFilterSeq seqFilter = new VariantFilterSeq(paramaters.getSeqFilterIn());
                varFilter.add(seqFilter);
            }
            if (paramaters.getMinMachR2() != 0.0d) {
                varFilter.add(new VariantFilterMachR2(paramaters.getMinMachR2()));
            }
            if(paramaters.filterOnAmbiguousSnps()){
                varFilter.add(new VariantFilterAmbigousSnp());
            }
        }

        SampleIdIncludeFilter sampleFilter = null;

        if (paramaters.getSampleFilterListFile() != null) {
            try {
                HashSet<String> samples = new HashSet<String>();

                BufferedReader sampleIdFilterReader = new BufferedReader(new FileReader(paramaters.getSampleFilterListFile()));
                String line;
                while ((line = sampleIdFilterReader.readLine()) != null) {
                    samples.add(line);
                }
                sampleFilter = new SampleIdIncludeFilter(samples);
            } catch (FileNotFoundException ex) {
                System.err.println("Unable to find file with samples to filter on at: " + paramaters.getSampleFilterListFile().getAbsolutePath());
                LOGGER.fatal("Unable to find file with samples to filter on at: " + paramaters.getSampleFilterListFile().getAbsolutePath());
                System.exit(1);
            } catch (IOException e) {
                System.err.println("Error reading file with samples to filter on at: " + paramaters.getSampleFilterListFile().getAbsolutePath() + " error: " + e.getMessage());
                LOGGER.fatal("Error reading file with samples to filter on at: " + paramaters.getSampleFilterListFile().getAbsolutePath(), e);
                System.exit(1);
            }
        }

        final RandomAccessGenotypeData inputData;

        try {
            inputData = paramaters.getInputType().createFilteredGenotypeData(paramaters.getInputBasePaths(), genotypeDataCache, varFilter, sampleFilter, paramaters.getForceSeqName(), paramaters.getMinimumPosteriorProbability());
        } catch (TabixFileNotFoundException e) {
            System.err.println("Tabix file not found for input data at: " + e.getPath() + "\n"
                    + "Please see README on how to create a tabix file");
            LOGGER.fatal("Tabix file not found for input data at: " + e.getPath() + "\n"
                    + "Please see README on how to create a tabix file");
            System.exit(1);
            return;
        } catch (IOException e) {
            System.err.println("Error reading input data: " + e.getMessage());
            LOGGER.fatal("Error reading input data: " + e.getMessage(), e);
            System.exit(1);
            return;
        } catch (IncompatibleMultiPartGenotypeDataException e) {
            System.err.println("Error combining the impute genotype data files: " + e.getMessage());
            LOGGER.fatal("Error combining the impute genotype data files: " + e.getMessage(), e);
            System.exit(1);
            return;
        } catch (GenotypeDataException e) {
            System.err.println("Error reading input data: " + e.getMessage());
            LOGGER.fatal("Error reading input data: " + e.getMessage(), e);
            System.exit(1);
            return;
        }

        if (inputData.getSamples().isEmpty()) {
            if (sampleFilter != null) {
                System.err.println("Error reading input data: No samples left after sample filter");
                LOGGER.fatal("Error reading input data: No samples left after sample filter");
            } else {
                System.err.println("Error reading input data: No samples are loaded");
                LOGGER.fatal("Error reading input data: No samples are loaded");
            }
            System.exit(1);
            return;
        }

        System.out.println(
                "Input data loaded");
        LOGGER.info(
                "Input data loaded");

        final RandomAccessGenotypeData refData;
        final ModifiableGenotypeData aligedInputData;
        if (paramaters.getRefBasePath() != null) {

            try {
                refData = paramaters.getRefType().createGenotypeData(paramaters.getRefBasePath(), genotypeDataCache);
            } catch (TabixFileNotFoundException e) {
                System.err.println("Tabix file not found for reference data at: " + e.getPath() + "\n"
                        + "Please see README on how to create a tabix file");
                LOGGER.fatal("Tabix file not found for reference data at: " + e.getPath() + "\n"
                        + "Please see README on how to create a tabix file");
                System.exit(1);
                return;
            } catch (IOException e) {
                System.err.println("Error reading reference data: " + e.getMessage());
                LOGGER.fatal("Error reading reference data: " + e.getMessage(), e);
                System.exit(1);
                return;
            } catch (IncompatibleMultiPartGenotypeDataException e) {
                System.err.println("Error combining the reference genotype data files: " + e.getMessage());
                LOGGER.fatal("Error combining the reference genotype data files: " + e.getMessage(), e);
                System.exit(1);
                return;
            } catch (GenotypeDataException e) {
                System.err.println("Error reading reference data: " + e.getMessage());
                LOGGER.fatal("Error reading reference data: " + e.getMessage(), e);
                System.exit(1);
                return;
            }

            System.out.println(
                    "Reference data loaded");
            LOGGER.info(
                    "Reference data loaded");

            Aligner aligner = new Aligner();

            try {
                System.out.println("Beginning alignment");
                aligedInputData = aligner.alignToRef(inputData, refData, paramaters.getMinLdToIncludeAlign(), paramaters.getMinSnpsToAlignOn(), paramaters.getFlankSnpsToConsider(), paramaters.isLdCheck(), paramaters.isUpdateId(), paramaters.isKeep(), paramaters.getSnpUpdateFile(), paramaters.getMaxMafForMafAlignment(), paramaters.getSnpLogFile(), paramaters.getMatchRefAllele());
            } catch (LdCalculatorException e) {
                System.err.println("Error in LD calculation: " + e.getMessage());
                LOGGER.fatal("Error in LD calculation: " + e.getMessage(), e);
                System.exit(1);
                return;
            } catch (GenotypeDataException e) {
                System.err.println("Error in alignment: " + e.getMessage());
                LOGGER.fatal("Error in alignment: " + e.getMessage(), e);
                System.exit(1);
                return;
            } catch (GenotypeAlignmentException e) {
                System.err.println("Error in alignment: " + e.getMessage());
                LOGGER.fatal("Error in alignment: " + e.getMessage(), e);
                System.exit(1);
                return;
            } catch (IOException e) {
                System.err.println("Error in alignment: " + e.getMessage());
                LOGGER.fatal("Error in alignment: " + e.getMessage(), e);
                System.exit(1);
                return;
            }

            System.out.println(
                    "Alignment complete");
            LOGGER.info(
                    "Alignment complete");

            System.out.println(
                    "Excluded in total " + DEFAULT_NUMBER_FORMATTER.format(aligedInputData.getExcludedVariantCount()) + " variants during alignment phase");
            LOGGER.info(
                    "Excluded in total " + DEFAULT_NUMBER_FORMATTER.format(aligedInputData.getExcludedVariantCount()) + " variants during alignment phase");
        } else {
            refData = null;
            aligedInputData = null;
            System.out.println("No reference specified. Do conversion without alignment");
            LOGGER.info("No reference specified. Do conversion without alignment");
        }

        if (paramaters.getInputType() == RandomAccessGenotypeDataReaderFormats.SHAPEIT2 && paramaters.getOutputType() == GenotypedDataWriterFormats.PLINK_BED) {
            System.out.println("WARNING: converting phased SHAPEIT2 data to binary Plink data. A BED file stores AB genotypes in the same manner as BA genotypes, thus all phasing will be lost.");
            LOGGER.warn("WARNING: converting phased SHAPEIT2 data to binary Plink data. A BED file stores AB genotypes in the same manner as BA genotypes, thus all phasing will be lost.");
        }

        if (paramaters.getOutputType() == GenotypedDataWriterFormats.GEN && !inputData.isOnlyContaingSaveProbabilityGenotypes()) {

            if (paramaters.getInputType() == RandomAccessGenotypeDataReaderFormats.VCF || paramaters.getInputType() == RandomAccessGenotypeDataReaderFormats.VCF_FOLDER) {
                System.out.println("WARNING: writing dosage genotype data to .gen posterior probabilities file. If sample does not have the GP field for a genotype then using heuristic method to convert to probabilities, this is not guaranteed to be accurate. See manual for more details.");
                LOGGER.warn("WARNING: writing dosage genotype data to .gen posterior probabilities file. If sample does not have the GP field for a genotype then using heuristic method to convert to probabilities, this is not guaranteed to be accurate. See manual for more details.");
            } else {
                System.out.println("WARNING: writing dosage genotype data to .gen posterior probabilities file. Using heuristic method to convert to probabilities, this is not guaranteed to be accurate. See manual for more details.");
                LOGGER.warn("WARNING: writing dosage genotype data to .gen posterior probabilities file. Using heuristic method to convert to probabilities, this is not guaranteed to be accurate. See manual for more details.");
            }

        }

        System.out.println(
                "Writing results");
        LOGGER.info(
                "Writing results");

        try {
            GenotypeWriter inputDataWriter = paramaters.getOutputType().createGenotypeWriter(aligedInputData == null ? inputData : aligedInputData);
            inputDataWriter.write(paramaters.getOutputBasePath());
        } catch (IOException e) {
            System.err.println("Error writing output data: " + e.getMessage());
            LOGGER.fatal("Error writing output data: " + e.getMessage(), e);
            System.exit(1);
            return;
        } catch (GenotypeDataException e) {
            System.err.println("Error writing output data: " + e.getMessage());
            LOGGER.fatal("Error writing output data: " + e.getMessage(), e);
            System.exit(1);
            return;
        }

        try {
            inputData.close();
            if (refData != null) {
                refData.close();
            }
        } catch (IOException ex) {
        }

        LOGGER.info(
                "Output data written");
        LOGGER.info(
                "Program complete");

        System.out.println(
                "Output data written");
        System.out.println(
                "Program complete");

    }
}
