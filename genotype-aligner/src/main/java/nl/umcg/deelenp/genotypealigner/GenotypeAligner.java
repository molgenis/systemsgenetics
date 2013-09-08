package nl.umcg.deelenp.genotypealigner;

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
import org.molgenis.genotype.GenotypeWriter;
import org.molgenis.genotype.GenotypedDataWriterFormats;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.modifiable.ModifiableGenotypeData;
import org.molgenis.genotype.multipart.IncompatibleMultiPartGenotypeDataException;
import org.molgenis.genotype.util.LdCalculatorException;

@SuppressWarnings("static-access")
class GenotypeAligner {

	private static final String VERSION = ResourceBundle.getBundle("verion").getString("application.version");
	private static final String HEADER =
			"  /---------------------------------------\\\n"
			+ "  |            Genotype Aligner           |\n"
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
	private static final int DEFAULT_FLANK_VARIANTS_TO_CONSIDER = 50;
	/**
	 * The lowest allowed minimum for the number of SNPs needed to align on
	 */
	private static final int MIN_MIN_VARIANTS_TO_ALIGN_ON = 3;
	/**
	 * The default minimum LD before using a SNP for LD alignment
	 */
	private static final double DEFAULT_MIN_LD_TO_INCLUDE_ALIGN = 0.3;

	static {

		LOGGER = Logger.getLogger(GenotypeAligner.class);

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
				.isRequired()
				.create("r");
		OPTIONS.addOption(option);

		option = OptionBuilder.withArgName("type")
				.hasArg()
				.withDescription("The input data type. \n"
				+ "* PED_MAP - plink PED MAP files.\n"
				+ "* PLINK_BED - plink BED BIM FAM files.\n"
				+ "* VCF - bgziped vcf with tabix index file\n"
				+ "* VCFFOLDER - matches all bgziped vcf files + tabix index in a folder\n"
				+ "* SHAPEIT2 - shapeit2 phased haplotypes .haps & .sample")
				.withLongOpt("inputType")
				.isRequired()
				.create("I");
		OPTIONS.addOption(option);

		option = OptionBuilder.withArgName("type")
				.hasArg()
				.withDescription("The reference data type. \n"
				+ "* PED_MAP - plink PED MAP files.\n"
				+ "* PLINK_BED - plink BED BIM FAM files.\n"
				+ "* VCF - bgziped vcf with tabix index file\n"
				+ "* VCF_FOLDER - matches all bgziped vcf files + tabix index in a folder\n"
				+ "* SHAPEIT2 - shapeit2 phased haplotypes .haps & .sample")
				.withLongOpt("refType")
				.isRequired()
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
				.withDescription("The output data type. \n"
				+ "* PED_MAP - plink PED MAP files. \n"
				+ "* PLINK_BED - plink BED BIM FAM files.\n"
				+ "* SHAPEIT2 - shapeit2 phased haplotypes.")
				.withLongOpt("outputType")
				.isRequired()
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
				.withDescription("Minum LD between variant to align or check and a flanking variants in both input as reference. Defaults to " + DEFAULT_MIN_LD_TO_INCLUDE_ALIGN + ". It is NOT recommend to set this to zero")
				.withLongOpt("min-ld")
				.create("l");
		OPTIONS.addOption(option);

		option = OptionBuilder.withArgName("int")
				.hasArg()
				.withDescription("Minimum number of variants above ld-cutoff to do LD alignment. Variants that do not meet this requerement are excluded. Defaults to " + DEFAULT_MIN_VARIANTS_TO_ALIGN_ON + ". Min value: " + MIN_MIN_VARIANTS_TO_ALIGN_ON)
				.withLongOpt("min-variants")
				.create("m");
		OPTIONS.addOption(option);

		option = OptionBuilder.withArgName("boolean")
				.withDescription("Check ld structure of all variants after alignment and exclude variants that deviate")
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
				.withDescription("Shapeit2 does not output the sequence name in the first column of the haplotype file. Use this option to force the chromosome for all variants. This option is only valid incombination with --inputType SHAPEIT2")
				.withLongOpt("forceChr")
				.create("f");
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
		System.out.println("More information: github.com/PatrickDeelen/GenotypeAligner/wiki");
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
			inputType = RandomAccessGenotypeDataReaderFormats.valueOf(commandLine.getOptionValue('I').toUpperCase());
		} catch (IllegalArgumentException e) {
			System.err.println("Error parsing --inputType \"" + commandLine.getOptionValue('I') + "\" is not a valid input data format");
			System.exit(1);
			return;
		}

		final String refBasePath = commandLine.getOptionValue('r');

		final RandomAccessGenotypeDataReaderFormats refType;
		try {
			refType = RandomAccessGenotypeDataReaderFormats.valueOf(commandLine.getOptionValue('R').toUpperCase());
		} catch (IllegalArgumentException e) {
			System.err.println("Error parsing --refType \"" + commandLine.getOptionValue('R') + "\" is not a valid reference data format");
			System.exit(1);
			return;
		}

		final String outputBasePath = commandLine.getOptionValue('o');

		final GenotypedDataWriterFormats outputType;
		try {
			outputType = GenotypedDataWriterFormats.valueOf(commandLine.getOptionValue('O').toUpperCase());
		} catch (IllegalArgumentException e) {
			System.err.println("Error parsing --outputType \"" + commandLine.getOptionValue('O') + "\" is not a valid output data format");
			System.exit(1);
			return;
		}

		final boolean debugMode = commandLine.hasOption('d');
		final boolean updateId = commandLine.hasOption("id");

		final int minSnpsToAlignOn;
		final int flankSnpsToConsider;
		final double minLdToIncludeAlign;

		try {
			minSnpsToAlignOn = commandLine.hasOption('v') ? Integer.parseInt(commandLine.getOptionValue('s')) : DEFAULT_MIN_VARIANTS_TO_ALIGN_ON;
		} catch (NumberFormatException e) {
			System.err.println("Error parsing --variants \"" + commandLine.getOptionValue('s') + "\" is not an int");
			System.exit(1);
			return;
		}

		try {
			flankSnpsToConsider = commandLine.hasOption('m') ? Integer.parseInt(commandLine.getOptionValue('m')) : DEFAULT_FLANK_VARIANTS_TO_CONSIDER;
		} catch (NumberFormatException e) {
			System.err.println("Error parsing --min-variants \"" + commandLine.getOptionValue('s') + "\" is not an int");
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
		
		final String forceSeqName;
		
	
		forceSeqName = commandLine.hasOption('f') ? commandLine.getOptionValue('f') : null;
		

		final boolean ldCheck = commandLine.hasOption('c');
		final boolean keep = commandLine.hasOption('k');

		File logFile = new File(outputBasePath + ".log");
		if (logFile.getParentFile() != null && !logFile.getParentFile().isDirectory()) {
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

		LOGGER.info("\n" + HEADER);
		LOGGER.info("Version: " + VERSION);
		LOGGER.info("Log level: " + LOGGER.getLevel());

		System.out.println("Started logging");
		System.out.println();

		printOptions(inputBasePath, inputType, refBasePath, refType, outputBasePath, outputType, minSnpsToAlignOn, flankSnpsToConsider, minLdToIncludeAlign, ldCheck, debugMode, updateId, keep, forceSeqName);


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
		
		if(forceSeqName != null && inputType != RandomAccessGenotypeDataReaderFormats.SHAPEIT2){
			System.err.println("Error cannot force sequence name of: " + inputType.getName());
			System.exit(1);
			return;
		}

		int genotypeDataCache = flankSnpsToConsider * 4;
		final RandomAccessGenotypeData inputData;

		try {
			inputData = inputType.createGenotypeData(inputBasePath, genotypeDataCache, forceSeqName);
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
		}

		System.out.println("Input data loaded");
		LOGGER.info("Input data loaded");

		final RandomAccessGenotypeData refData;
		try {
			refData = refType.createGenotypeData(refBasePath, genotypeDataCache);
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
		}
		System.out.println("Reference data loaded");
		LOGGER.info("Reference data loaded");

		Aligner aligner = new Aligner();

		ModifiableGenotypeData aligedInputData;

		try {
			System.out.println("Beginning alignment");
			aligedInputData = aligner.alignToRef(inputData, refData, minLdToIncludeAlign, minSnpsToAlignOn, flankSnpsToConsider, ldCheck, updateId, keep);
		} catch (LdCalculatorException e) {
			System.err.println("Error in LD caculation" + e.getMessage());
			LOGGER.fatal("Error in LD caculation" + e.getMessage(), e);
			System.exit(1);
			return;
		}

		System.out.println("Alignment complete");
		LOGGER.info("Alignment complete");

		System.out.println("Excluded in total " + aligedInputData.getExcludedVariantCount() + " variants");
		LOGGER.info("Excluded in total " + aligedInputData.getExcludedVariantCount() + " variants");

		System.out.println("Writing results");
		LOGGER.info("Writing results");


		GenotypeWriter inputDataWriter = outputType.createGenotypeWriter(aligedInputData);

		try {
			inputDataWriter.write(outputBasePath);
		} catch (IOException e) {
			System.err.println("Error writing output data: " + e.getMessage());
			LOGGER.fatal("Error writing output data: " + e.getMessage(), e);
			System.exit(1);
			return;
		}
		try {
			inputData.close();
			refData.close();
		} catch (IOException ex) {
			
		}
				
		LOGGER.info("Output data writen");
		LOGGER.info("Program complete");

		System.out.println("Output data writen");
		System.out.println("Program complete");

	}

	private static void printOptions(String inputBasePath, RandomAccessGenotypeDataReaderFormats inputType, String refBasePath, RandomAccessGenotypeDataReaderFormats refType, String outputBasePath, GenotypedDataWriterFormats outputType, int minSnpsToAlignOn, int flankSnpsToConsider, double minLdToIncludeAlign, boolean ldCheck, boolean debugMode, boolean updateId, boolean keep, String forceSeqName) {


		System.out.println("Interpreted arugments: ");
		System.out.println(" - Input base path: " + inputBasePath);
		LOGGER.info("Input base path: " + inputBasePath);
		System.out.println(" - Input data type: " + inputType.getName());
		LOGGER.info("Input data type: " + inputType.getName());
		System.out.println(" - Reference base path: " + refBasePath);
		LOGGER.info("Reference base path: " + refBasePath);
		System.out.println(" - Reference data type: " + refType.getName());
		LOGGER.info("Reference data type: " + refType.getName());
		System.out.println(" - Output base path: " + outputBasePath);
		LOGGER.info("Output base path: " + outputBasePath);
		System.out.println(" - Output data type: " + outputType.getName());
		LOGGER.info("Output data type: " + outputType.getName());
		System.out.println(" - Number of flank variants to consider for LD alignment: " + flankSnpsToConsider);
		LOGGER.info("Number of flank variants to consider for LD alignment: " + flankSnpsToConsider);
		System.out.println(" - Minimum LD of flanking variants before using for LD alignment: " + minLdToIncludeAlign);
		LOGGER.info("Minimum LD of flanking variants before using for LD alignment: " + minLdToIncludeAlign);
		System.out.println(" - Minimum number of variants needed to for LD aligment: " + minSnpsToAlignOn);
		LOGGER.info("Minimum number of variants needed to for LD aligment: " + minSnpsToAlignOn);
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
