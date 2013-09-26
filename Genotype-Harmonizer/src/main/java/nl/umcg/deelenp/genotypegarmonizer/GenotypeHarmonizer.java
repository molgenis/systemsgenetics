package nl.umcg.deelenp.genotypegarmonizer;

import java.io.File;
import java.io.IOException;
import java.util.ResourceBundle;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
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
import org.molgenis.genotype.tabix.TabixFileNotFoundException;
import org.molgenis.genotype.util.LdCalculatorException;

/**
 * 
 * @author Patrick Deelen
 */
@SuppressWarnings("static-access")
class GenotypeHarmonizer {

	private static final String VERSION = ResourceBundle.getBundle("verion").getString("application.version");
	private static final String HEADER =
			"  /---------------------------------------\\\n"
			+ "  |          Genotype Harmonizer          |\n"
			+ "  |                                       |\n"
			+ "  |             Patrick Deelen            |\n"
			+ "  |        patrickdeelen@gmail.com        |\n"
			+ "  |                                       |\n"
			+ "  | Harm-Jan Westra, Joeri van der Velde, |\n"
			+ "  |    Marc Jan Bonder, Erwin Winder,     |\n"
			+ "  |      Lude Franke, Morris Swertz       |\n"
			+ "  |                                       |\n"
			+ "  |     Genomics Coordication Center      |\n"
			+ "  |  University Medical Center Groningen  |\n"
			+ "  \\---------------------------------------/";
	private static final Logger LOGGER;
	private static final Options OPTIONS;
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
	 * The lowest allowed minimum for the number of SNPs needed to align on
	 */
	private static final int MIN_MIN_VARIANTS_TO_ALIGN_ON = 3;
	/**
	 * The default minimum LD before using a SNP for LD alignment
	 */
	private static final double DEFAULT_MIN_LD_TO_INCLUDE_ALIGN = 0.3;
	/**
	 * The default maximum LD before the minor allele frequency (MAF) is used as backup for alignment
	 */
	private static final double DEFAULT_MAX_MAF_FOR_MAF_ALIGNMENT = 0;

	static {

		LOGGER = Logger.getLogger(GenotypeHarmonizer.class);

		OPTIONS = new Options();

		Option option;

		option = OptionBuilder.withArgName("basePath")
				.hasArg()
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
				+ "* SHAPEIT2 - shapeit2 phased haplotypes .haps & .sample")
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
				+ "* SHAPEIT2 - shapeit2 phased haplotypes .haps & .sample")
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
				+ "* PED_MAP - plink PED MAP files. \n"
				+ "* PLINK_BED - plink BED BIM FAM files.\n"
				+ "* SHAPEIT2 - shapeit2 phased haplotypes.")
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
				.withDescription("Minimum number of variants above ld-cutoff to do LD alignment. Variants that do not meet this requirement are excluded. Defaults to " + DEFAULT_MIN_VARIANTS_TO_ALIGN_ON + ". Min value: " + MIN_MIN_VARIANTS_TO_ALIGN_ON)
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

	}

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

		System.out.flush(); //flush to make sure header is before errors
		Thread.sleep(25); //Allows flush to complete

		if (args.length == 0) {
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(" ", OPTIONS);
			System.exit(1);
		}

		final CommandLineParser parser = new PosixParser();
		final CommandLine commandLine;
		try {
			commandLine = parser.parse(OPTIONS, args, true);
		} catch (ParseException ex) {
			System.err.println("Invalid command line arguments: ");
			System.err.println(ex.getMessage());
			System.err.println();
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(" ", OPTIONS);
			System.exit(1);
			return;
		}


		final String inputBasePath = commandLine.getOptionValue('i');
		final RandomAccessGenotypeDataReaderFormats inputType;

		try {
			if (commandLine.hasOption('I')) {
				inputType = RandomAccessGenotypeDataReaderFormats.valueOf(commandLine.getOptionValue('I').toUpperCase());
			} else {
				if (inputBasePath.endsWith(".vcf")) {
					System.err.println("Only vcf.gz is supported. Please see manual on how to do create a vcf.gz file.");
				}
				try {
					inputType = RandomAccessGenotypeDataReaderFormats.matchFormatToPath(inputBasePath);
				} catch (GenotypeDataException e) {
					System.err.println("Unable to determine input type based on specified path. Please specify --inputType");
					System.exit(1);
					return;
				}
			}
		} catch (IllegalArgumentException e) {
			System.err.println("Error parsing --inputType \"" + commandLine.getOptionValue('I') + "\" is not a valid input data format");
			System.exit(1);
			return;
		}

		final String refBasePath;

		final RandomAccessGenotypeDataReaderFormats refType;

		if (commandLine.hasOption('r')) {

			refBasePath = commandLine.getOptionValue('r');

			try {
				if (commandLine.hasOption('R')) {
					refType = RandomAccessGenotypeDataReaderFormats.valueOf(commandLine.getOptionValue('R').toUpperCase());
				} else {
					if (refBasePath.endsWith(".vcf")) {
						System.err.println("Only vcf.gz is supported. Please see manual on how to do create a vcf.gz file.");
						System.exit(1);
						return;
					}
					try {
						refType = RandomAccessGenotypeDataReaderFormats.matchFormatToPath(refBasePath);
					} catch (GenotypeDataException e) {
						System.err.println("Unable to determine reference type based on specified path. Please specify --refType");
						System.exit(1);
						return;
					}
				}

			} catch (IllegalArgumentException e) {
				System.err.println("Error parsing --refType \"" + commandLine.getOptionValue('R') + "\" is not a valid reference data format");
				System.exit(1);
				return;
			}

		} else {
			refBasePath = null;
			refType = null;
		}




		final String outputBasePath = commandLine.getOptionValue('o');

		GenotypedDataWriterFormats outputType;

		if (commandLine.hasOption('O')) {
			try {
				outputType = GenotypedDataWriterFormats.valueOf(commandLine.getOptionValue('O').toUpperCase());
			} catch (IllegalArgumentException e) {
				System.err.println("Error parsing --outputType \"" + commandLine.getOptionValue('O') + "\" is not a valid output data format");
				System.exit(1);
				return;
			}
		} else {
			try {
				outputType = GenotypedDataWriterFormats.valueOf(inputType.name());
			} catch (IllegalArgumentException e) {
				outputType = GenotypedDataWriterFormats.PLINK_BED;
			}
		}


		final boolean debugMode = commandLine.hasOption('d');
		final boolean updateId = commandLine.hasOption("id");
		final int minSnpsToAlignOn;
		final int flankSnpsToConsider;
		final double minLdToIncludeAlign;
		final double maxMafForMafAlignment;


		try {
			minSnpsToAlignOn = commandLine.hasOption('m') ? Integer.parseInt(commandLine.getOptionValue('m')) : DEFAULT_MIN_VARIANTS_TO_ALIGN_ON;
		} catch (NumberFormatException e) {
			System.err.println("Error parsing --min-variants \"" + commandLine.getOptionValue('m') + "\" is not an int");
			System.exit(1);
			return;
		}

		try {
			flankSnpsToConsider = commandLine.hasOption('v') ? Integer.parseInt(commandLine.getOptionValue('v')) : DEFAULT_FLANK_VARIANTS_TO_CONSIDER;
		} catch (NumberFormatException e) {
			System.err.println("Error parsing --variants \"" + commandLine.getOptionValue('v') + "\" is not an int");
			System.exit(1);
			return;
		}

		try {
			minLdToIncludeAlign = commandLine.hasOption('l') ? Double.parseDouble(commandLine.getOptionValue('l')) : DEFAULT_MIN_LD_TO_INCLUDE_ALIGN;
		} catch (NumberFormatException e) {
			System.err.println("Error parsing --min-ld \"" + commandLine.getOptionValue('s') + "\" is not a double");
			System.exit(1);
			return;
		}
		
		try{
			maxMafForMafAlignment = commandLine.hasOption("ma") ? Double.parseDouble(commandLine.getOptionValue("ma")) : DEFAULT_MAX_MAF_FOR_MAF_ALIGNMENT;
		} catch (NumberFormatException e) {
			System.err.println("Error parsing --mafAlign \"" + commandLine.getOptionValue("ma") + "\" is not a double");
			System.exit(1);
			return;
		}
		
		
		final String forceSeqName = commandLine.hasOption('f') ? commandLine.getOptionValue('f') : null;
		final boolean ldCheck = commandLine.hasOption('c');
		final boolean keep = commandLine.hasOption('k');
		File logFile = new File(outputBasePath + ".log");
		File snpUpdateFile = new File(outputBasePath + ".idUpdates");

		if (logFile.getParentFile()
				!= null && !logFile.getParentFile().isDirectory()) {
			if (!logFile.getParentFile().mkdirs()) {
				System.err.println("Failed to create output folder: " + logFile.getParent());
				System.exit(1);
			}
		}


		try {
			FileAppender logAppender = new FileAppender(new SimpleLayout(), logFile.getCanonicalPath(), false);
			LOGGER.getRootLogger().removeAllAppenders();
			LOGGER.getRootLogger().addAppender(logAppender);
			if (debugMode) {
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
		LOGGER.info(
				"Version: " + VERSION);
		LOGGER.info(
				"Log level: " + LOGGER.getLevel());

		System.out.println(
				"Started logging");
		System.out.println();

		printOptions(inputBasePath, inputType, refBasePath, refType, outputBasePath, outputType, minSnpsToAlignOn, flankSnpsToConsider, minLdToIncludeAlign, ldCheck, debugMode, updateId, keep, forceSeqName, maxMafForMafAlignment);
		
		if (minSnpsToAlignOn < MIN_MIN_VARIANTS_TO_ALIGN_ON) {
			LOGGER.fatal("the specified --min-variants < " + MIN_MIN_VARIANTS_TO_ALIGN_ON);
			System.err.println("the specified --min-variants < " + MIN_MIN_VARIANTS_TO_ALIGN_ON);
			System.exit(1);
			return;
		}
		if (flankSnpsToConsider < minLdToIncludeAlign) {
			LOGGER.fatal("--variants < --min-variants");
			System.err.println("--variants < --min-variants");
			System.exit(1);
			return;
		}

		if (inputBasePath.equals(refBasePath)) {
			LOGGER.fatal("Study data and reference data cannot be the same data");
			System.err.println("Study data and reference data cannot be the same data");
			System.exit(1);
			return;
		}

		if (inputBasePath.equals(outputBasePath)) {
			LOGGER.fatal("Study input can not be the same as output");
			System.err.println("Study input can not be the same as output");
			System.exit(1);
			return;
		}
		if (forceSeqName != null && inputType != RandomAccessGenotypeDataReaderFormats.SHAPEIT2) {
			System.err.println("Error cannot force sequence name of: " + inputType.getName());
			System.exit(1);
			return;
		}
		int genotypeDataCache = flankSnpsToConsider * 4;
		final RandomAccessGenotypeData inputData;


		try {
			inputData = inputType.createGenotypeData(inputBasePath, genotypeDataCache, forceSeqName);
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

		System.out.println(
				"Input data loaded");
		LOGGER.info(
				"Input data loaded");


		final RandomAccessGenotypeData refData;
		final ModifiableGenotypeData aligedInputData;
		if (refBasePath != null) {

			try {
				refData = refType.createGenotypeData(refBasePath, genotypeDataCache);
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
			if (inputType == RandomAccessGenotypeDataReaderFormats.SHAPEIT2 && outputType == GenotypedDataWriterFormats.PLINK_BED) {
				System.out.println("WARNING: converting phased SHAPEIT2 data to binary Plink data. A BED file stores AB genotypes in the same manner as BA genotypes, thus all phasing will be lost.");
				LOGGER.warn("WARNING: converting phased SHAPEIT2 data to binary Plink data. A BED file stores AB genotypes in the same manner as BA genotypes, thus all phasing will be lost.");
			}


			try {
				System.out.println("Beginning alignment");
				aligedInputData = aligner.alignToRef(inputData, refData, minLdToIncludeAlign, minSnpsToAlignOn, flankSnpsToConsider, ldCheck, updateId, keep, snpUpdateFile, maxMafForMafAlignment);
			} catch (LdCalculatorException e) {
				System.err.println("Error in LD calculation" + e.getMessage());
				LOGGER.fatal("Error in LD calculation" + e.getMessage(), e);
				System.exit(1);
				return;
			} catch (GenotypeDataException e) {
				System.err.println("Error in alignment" + e.getMessage());
				LOGGER.fatal("Error in alignment" + e.getMessage(), e);
				System.exit(1);
				return;
			} catch (IOException e) {
				System.err.println("Error in alignment" + e.getMessage());
				LOGGER.fatal("Error in alignment" + e.getMessage(), e);
				System.exit(1);
				return;
			}

			System.out.println(
					"Alignment complete");
			LOGGER.info(
					"Alignment complete");

			System.out.println(
					"Excluded in total " + aligedInputData.getExcludedVariantCount() + " variants");
			LOGGER.info(
					"Excluded in total " + aligedInputData.getExcludedVariantCount() + " variants");
		} else {
			refData = null;
			aligedInputData = null;
			System.out.println("No reference specified. Do conversion without alignment");
			LOGGER.info("No reference specified. Do conversion without alignment");		
		}

		System.out.println(
				"Writing results");
		LOGGER.info(
				"Writing results");


		try {
			GenotypeWriter inputDataWriter = outputType.createGenotypeWriter(aligedInputData == null ? inputData : aligedInputData);
			inputDataWriter.write(outputBasePath);
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

	private static void printOptions(String inputBasePath, RandomAccessGenotypeDataReaderFormats inputType, String refBasePath, RandomAccessGenotypeDataReaderFormats refType, String outputBasePath, GenotypedDataWriterFormats outputType, int minSnpsToAlignOn, int flankSnpsToConsider, double minLdToIncludeAlign, boolean ldCheck, boolean debugMode, boolean updateId, boolean keep, String forceSeqName, double maxMafForMafAlignment) {


		System.out.println("Interpreted arguments: ");
		System.out.println(" - Input base path: " + inputBasePath);
		LOGGER.info("Input base path: " + inputBasePath);
		System.out.println(" - Input data type: " + inputType.getName());
		LOGGER.info("Input data type: " + inputType.getName());
		System.out.println(" - Reference base path: " + (refBasePath == null ? "no reference set" : refBasePath));
		LOGGER.info("Reference base path: " + (refBasePath == null ? "no reference set" : refBasePath));
		System.out.println(" - Reference data type: " + (refBasePath == null ? "no reference set" : refType.getName()));
		LOGGER.info("Reference data type: " + (refBasePath == null ? "no reference set" : refType.getName()));
		System.out.println(" - Output base path: " + outputBasePath);
		LOGGER.info("Output base path: " + outputBasePath);
		System.out.println(" - Output data type: " + outputType.getName());
		LOGGER.info("Output data type: " + outputType.getName());
		System.out.println(" - Number of flank variants to consider for LD alignment: " + flankSnpsToConsider);
		LOGGER.info("Number of flank variants to consider for LD alignment: " + flankSnpsToConsider);
		System.out.println(" - Minimum LD of flanking variants before using for LD alignment: " + minLdToIncludeAlign);
		LOGGER.info("Minimum LD of flanking variants before using for LD alignment: " + minLdToIncludeAlign);
		System.out.println(" - Minimum number of variants needed to for LD alignment: " + minSnpsToAlignOn);
		LOGGER.info("Minimum number of variants needed to for LD alignment: " + minSnpsToAlignOn);
		System.out.println(" - Maximum MAF of variants to use minor allele as backup for alignment: " + maxMafForMafAlignment);
		LOGGER.info("Maximum MAF of variants to use minor allele as backup for alignment: " + maxMafForMafAlignment);
		System.out.println(" - LD checker " + (ldCheck ? "on" : "off"));
		LOGGER.info("LD checker " + (ldCheck ? "on" : "off"));
		System.out.println(" - Update study IDs: " + (updateId ? "yes" : "no"));
		LOGGER.info("Update study variant IDs: " + (updateId ? "yes" : "no"));
		LOGGER.info("Keep variants not in reference data: " + (keep ? "yes" : "no"));
		System.out.println(" - Keep variants not in reference data: " + (keep ? "yes" : "no"));
		LOGGER.info("Force input sequence name: " + (forceSeqName == null ? "not forcing" : "forcing to: " + forceSeqName));
		System.out.println(" - Force input sequence name: " + (forceSeqName == null ? "not forcing" : "forcing to: " + forceSeqName));
		LOGGER.info("Debug mode: " + (debugMode ? "on" : "off"));


		System.out.println();


	}
}
