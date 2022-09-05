/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.toolbox;

import htsjdk.tribble.SimpleFeature;
import java.io.File;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.log4j.Logger;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;

/**
 *
 * @author patri
 */
public class Options {

	private static final org.apache.commons.cli.Options OPTIONS;
	private static final Logger LOGGER = Logger.getLogger(Options.class);
	private static final SimpleFeature HLA = new SimpleFeature("6", 20000000, 40000000);

	private final Mode mode;

	private final String inputPath;
	private final File outputBasePath;
	private final File logFile;
	private final boolean debugMode;

	private final String[] genotypeBasePath;
	private final RandomAccessGenotypeDataReaderFormats genotypeType;
	private final File genotypeSamplesFile;
	private final File variantFilterFile;
	private final double mafFilter;

	private final File gwasCatalogFile;

	static {

		OPTIONS = new org.apache.commons.cli.Options();

		OptionBuilder.withArgName("string");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("On of the following modes:\n"
				+ "* OVERLAP_GWAS"
		);
		OptionBuilder.withLongOpt("mode");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("m"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Genetic primary input path for differnt modes");
		OptionBuilder.withLongOpt("input");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("i"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("The output path");
		OptionBuilder.withLongOpt("output");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("o"));

		OptionBuilder.withArgName("boolean");
		OptionBuilder.withDescription("Activate debug mode.");
		OptionBuilder.withLongOpt("debug");
		OPTIONS.addOption(OptionBuilder.create("d"));

		OptionBuilder.withArgName("basePath");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("The reference genotypes");
		OptionBuilder.withLongOpt("referenceGenotypes");
		OPTIONS.addOption(OptionBuilder.create("r"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Samples to include from reference genotypes");
		OptionBuilder.withLongOpt("referenceSamples");
		OPTIONS.addOption(OptionBuilder.create("rs"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("File with variants to include in gene p-value calculation (optional)");
		OptionBuilder.withLongOpt("variantFilter");
		OPTIONS.addOption(OptionBuilder.create("vf"));

		OptionBuilder.withArgName("type");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("The reference genotype data type. If not defined will attempt to automatically select the first matching dataset on the specified path\n"
				+ "* PED_MAP - plink PED MAP files.\n"
				+ "* PLINK_BED - plink BED BIM FAM files.\n"
				+ "* VCF - bgziped vcf with tabix index file\n"
				+ "* VCFFOLDER - matches all bgziped vcf files + tabix index in a folder\n"
				+ "* SHAPEIT2 - shapeit2 phased haplotypes .haps & .sample\n"
				+ "* GEN - Oxford .gen & .sample\n"
				+ "* BGEN - Oxford .bgen & optionally .sample\n"
				+ "* TRITYPER - TriTyper format folder");
		OptionBuilder.withLongOpt("referenceGenotypeFormat");
		OPTIONS.addOption(OptionBuilder.create("R"));
		
		OptionBuilder.withArgName("double");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Minimum MAF");
		OptionBuilder.withLongOpt("maf");
		OPTIONS.addOption(OptionBuilder.create("maf"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("GWAS catalog file");
		OptionBuilder.withLongOpt("gwasCatalog");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("gc"));

	}

	public Options(String... args) throws ParseException {

		final CommandLineParser parser = new PosixParser();
		final CommandLine commandLine = parser.parse(OPTIONS, args, false);

		inputPath = commandLine.getOptionValue('i');
		outputBasePath = new File(commandLine.getOptionValue('o'));
		logFile = new File(outputBasePath + ".log");
		debugMode = commandLine.hasOption('d');

		try {

			String modeString = commandLine.getOptionValue("m").toUpperCase();
			mode = Mode.valueOf(modeString);

		} catch (IllegalArgumentException e) {
			throw new ParseException("Error parsing --mode \"" + commandLine.getOptionValue("m") + "\" is not a valid mode");
		}

		if (mode.isNeedGenotypes()) {
			if (!commandLine.hasOption('r')) {

				throw new ParseException("--referenceGenotypes not specified");

			} else {

				genotypeBasePath = commandLine.getOptionValues('r');

				if (commandLine.hasOption("rs")) {
					genotypeSamplesFile = new File(commandLine.getOptionValue("rs"));
				} else {
					genotypeSamplesFile = null;
				}

				if (commandLine.hasOption("vf")) {
					variantFilterFile = new File(commandLine.getOptionValue("vf"));
				} else {
					variantFilterFile = null;
				}
				
				if (commandLine.hasOption("maf")) {
					try {
						mafFilter = Double.parseDouble(commandLine.getOptionValue("maf"));
					} catch (NumberFormatException e) {
						throw new ParseException("Error parsing --maf \"" + commandLine.getOptionValue("maf") + "\" is not an double");
					}
				} else {
					mafFilter = 0;
				}

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
							throw new ParseException("Unable to determine reference type based on specified path. Please specify --refType");
						}
					}

				} catch (IllegalArgumentException e) {
					throw new ParseException("Error parsing --refType \"" + commandLine.getOptionValue('R') + "\" is not a valid reference data format");
				}

			}
		} else {
			genotypeBasePath = null;
			genotypeType = null;
			genotypeSamplesFile = null;
			variantFilterFile = null;
			mafFilter = 0;
		}

		if (mode.isNeedGwasCatalog()) {
			if (!commandLine.hasOption("gc")) {
				throw new ParseException("Please provide --gwasCatalog for mode: " + mode.name());
			}

			gwasCatalogFile = new File(commandLine.getOptionValue("gc"));

		} else {
			gwasCatalogFile = null;
		}

	}

	public static void printHelp() {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp(" ", OPTIONS);
	}

	public void printOptions() {

		LOGGER.info("Supplied options:");
		LOGGER.info(" * Mode: " + mode.name());
		LOGGER.info(" * Input path: " + inputPath);
		LOGGER.info(" * Ouput path: " + outputBasePath.getAbsolutePath());
		LOGGER.info(" * Debug mode: " + (debugMode ? "on" : "off"));

		if (mode.isNeedGenotypes()) {

			StringBuilder genotypeBasePaths = new StringBuilder();
			for (String path : genotypeBasePath) {
				genotypeBasePaths.append(new File(path).getAbsolutePath());
				genotypeBasePaths.append(' ');
			}
			LOGGER.info(" * Reference genotype data: " + genotypeBasePaths);
			LOGGER.info(" * Reference genotype data type: " + genotypeType.getName());

			if (variantFilterFile != null) {
				LOGGER.info(" * Confining genotype data to variants in this file: " + variantFilterFile.getAbsolutePath());
			}
			if (genotypeSamplesFile != null) {
				LOGGER.info(" * Confining genotype data to samples in this file: " + genotypeSamplesFile.getAbsolutePath());
			}
		}

		if(mode.isNeedGwasCatalog()){
			LOGGER.info(" * GWAS catalog file: " + gwasCatalogFile.getAbsolutePath());
		}
		
	}

	public String[] getGenotypeBasePath() {
		return genotypeBasePath;
	}

	public RandomAccessGenotypeDataReaderFormats getGenotypeType() {
		return genotypeType;
	}

	public String getInputPath() {
		return inputPath;
	}

	public String getOutputBasePath() {
		return outputBasePath.getPath();
	}

	public File getLogFile() {
		return logFile;
	}

	public File getGenotypeSamplesFile() {
		return genotypeSamplesFile;
	}

	public File getVariantFilterFile() {
		return variantFilterFile;
	}

	public boolean isDebugMode() {
		return debugMode;
	}

	public Mode getMode() {
		return mode;
	}

	public File getGwasCatalogFile() {
		return gwasCatalogFile;
	}

	public double getMafFilter() {
		return mafFilter;
	}
	
}
