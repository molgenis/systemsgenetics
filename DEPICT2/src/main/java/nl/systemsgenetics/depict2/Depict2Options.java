/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import java.io.File;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.log4j.Logger;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;

/**
 *
 * @author patri
 */
public class Depict2Options {

	private static final Options OPTIONS;
	private static int numberOfThreadsToUse = Runtime.getRuntime().availableProcessors();
	private static final Logger LOGGER = Logger.getLogger(Depict2Options.class);

	private final Depict2Mode mode;

	private final String[] genotypeBasePath;
	private final RandomAccessGenotypeDataReaderFormats genotypeType;
	private final File genotypeSamplesFile;
	private final String outputBasePath;
	private final File geneInfoFile;
	private final String gwasZscoreMatrixPath;
	private final int numberOfPermutations;
	private final int windowExtend;
	private final double maxRBetweenVariants;
	private final File logFile;
	private final boolean debugMode;

	public boolean isDebugMode() {
		return debugMode;
	}

	static {

		OPTIONS = new Options();

		OptionBuilder.withArgName("string");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("On of the following modes:\n"
				+ "* RUN - Run the DEPICT2 prioritization.\n"
				+ "* CONVERT_TXT - Convert a txt z-score matrix to binary. Use --gwas and --output\n"
				+ "* CONVERT_EQTL - Convert binary matrix with eQTL z-scores from our pipeline. Use --gwas and --output");
		OptionBuilder.withLongOpt("mode");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("m"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("GWAS Z-sccore binary matrix. Rows variants, Cols phenotypes. Without .dat");
		OptionBuilder.withLongOpt("gwas");
		OPTIONS.addOption(OptionBuilder.create("g"));

		OptionBuilder.withArgName("basePath");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("The reference genotypes");
		OptionBuilder.withLongOpt("referenceGenotypes");
		OPTIONS.addOption(OptionBuilder.create("r"));

		OptionBuilder.withArgName("basePath");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Samples to include from reference genotypes");
		OptionBuilder.withLongOpt("referenceSamples");
		OPTIONS.addOption(OptionBuilder.create("rs"));

		OptionBuilder.withArgName("type");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("The reference genotype data type. If not defined will attempt to automatically select the first matching dataset on the specified path\n"
				+ "* PED_MAP - plink PED MAP files.\n"
				+ "* PLINK_BED - plink BED BIM FAM files.\n"
				+ "* VCF - bgziped vcf with tabix index file\n"
				+ "* VCFFOLDER - matches all bgziped vcf files + tabix index in a folder\n"
				+ "* SHAPEIT2 - shapeit2 phased haplotypes .haps & .sample\n"
				+ "* GEN - Oxford .gen & .sample\n"
				+ "* TRITYPER - TriTyper format folder");
		OptionBuilder.withLongOpt("referenceGenotypeFormat");
		OPTIONS.addOption(OptionBuilder.create("R"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("The output path");
		OptionBuilder.withLongOpt("output");
		OPTIONS.addOption(OptionBuilder.create("o"));

		OptionBuilder.withArgName("int");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Maximum number of calculation threads");
		OptionBuilder.withLongOpt("threads");
		OPTIONS.addOption(OptionBuilder.create("t"));

		OptionBuilder.withArgName("int");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Number of permutations");
		OptionBuilder.withLongOpt("permutations");
		OPTIONS.addOption(OptionBuilder.create("p"));

		OptionBuilder.withArgName("int");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Number of bases to add left and right of gene window");
		OptionBuilder.withLongOpt("window");
		OPTIONS.addOption(OptionBuilder.create("w"));

		OptionBuilder.withArgName("double");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Max correlation between variants to use (recommend = 0.95)");
		OptionBuilder.withLongOpt("variantCorrelation");
		OPTIONS.addOption(OptionBuilder.create("v"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("File with gene info. col1: geneName (ensg) col2: chr col3: startPos col4: stopPos");
		OptionBuilder.withLongOpt("genes");
		OPTIONS.addOption(OptionBuilder.create("ge"));

		OptionBuilder.withArgName("boolean");
		OptionBuilder.withDescription("Activate debug mode. This will result in a more verbose log file");
		OptionBuilder.withLongOpt("debug");
		OPTIONS.addOption(OptionBuilder.create("d"));

	}

	public Depict2Options(String... args) throws ParseException {

		final CommandLineParser parser = new PosixParser();
		final CommandLine commandLine = parser.parse(OPTIONS, args, false);

		if (commandLine.hasOption('t')) {
			try {
				numberOfThreadsToUse = Integer.parseInt(commandLine.getOptionValue('t'));
			} catch (NumberFormatException e) {
				throw new ParseException("Error parsing --threads \"" + commandLine.getOptionValue('t') + "\" is not an int");
			}
		}

		outputBasePath = commandLine.getOptionValue('o');
		logFile = new File(outputBasePath + ".log");
		debugMode = commandLine.hasOption('d');

		try {
			mode = Depict2Mode.valueOf(commandLine.getOptionValue("m"));
		} catch (IllegalArgumentException e) {
			throw new ParseException("Error parsing --mode \"" + commandLine.getOptionValue("m") + "\" is not a valid mode");
		}

		if (mode == Depict2Mode.CONVERT_TXT || mode == Depict2Mode.RUN || mode == Depict2Mode.CONVERT_EQTL) {

			if (!commandLine.hasOption("g")) {
				throw new ParseException("Please provide --gwas for mode: " + mode.name());
			}

			gwasZscoreMatrixPath = commandLine.getOptionValue('g');
		} else {
			gwasZscoreMatrixPath = null;
		}

		if (mode == Depict2Mode.RUN) {

			if (!commandLine.hasOption("p")) {
				throw new ParseException("--permutations not specified");
			} else {

				try {
					numberOfPermutations = Integer.parseInt(commandLine.getOptionValue('p'));
				} catch (NumberFormatException e) {
					throw new ParseException("Error parsing --permutations \"" + commandLine.getOptionValue('p') + "\" is not an int");
				}
			}

			if (!commandLine.hasOption("w")) {
				throw new ParseException("--window not specified");
			} else {
				try {
					windowExtend = Integer.parseInt(commandLine.getOptionValue('w'));
				} catch (NumberFormatException e) {
					throw new ParseException("Error parsing --window \"" + commandLine.getOptionValue('w') + "\" is not an int");
				}
			}

			if (!commandLine.hasOption("v")) {
				throw new ParseException("--variantCorrelation not specified");
			} else {
				try {
					maxRBetweenVariants = Double.parseDouble(commandLine.getOptionValue('v'));
				} catch (NumberFormatException e) {
					throw new ParseException("Error parsing --variantCorrelation \"" + commandLine.getOptionValue('v') + "\" is not an double");
				}
			}

			if (!commandLine.hasOption("ge")) {
				throw new ParseException("--genes not specified");
			} else {
				geneInfoFile = new File(commandLine.getOptionValue("ge"));
			}

			if (!commandLine.hasOption('r')) {

				throw new ParseException("--referenceGenotypes not specified");

			} else {

				genotypeBasePath = commandLine.getOptionValues('r');

				if (commandLine.hasOption("rs")) {
					genotypeSamplesFile = new File(commandLine.getOptionValue("rs"));
				} else {
					genotypeSamplesFile = null;
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
			geneInfoFile = null;
			maxRBetweenVariants = 0d;
			numberOfPermutations = 0;
			windowExtend = 0;
		}

	}

	public static void printHelp() {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp(" ", OPTIONS);
	}

	public void printOptions() {

		LOGGER.info("Supplied options:");
		
		LOGGER.info(" * Mode: " + mode.name());
		LOGGER.info(" * Ouput path: " + outputBasePath);

		switch (mode) {
			case CONVERT_EQTL:
				LOGGER.info(" * eQTL Z-score matrix: " + gwasZscoreMatrixPath);
				break;
			case CONVERT_TXT:
				LOGGER.info(" * Gwas Z-score matrix: " + gwasZscoreMatrixPath);
				break;
			case RUN:
				LOGGER.info(" * Gwas Z-score matrix: " + gwasZscoreMatrixPath);

				if (genotypeBasePath != null) {
					StringBuilder genotypeBasePaths = new StringBuilder();
					for (String path : genotypeBasePath) {
						genotypeBasePaths.append(path);
						genotypeBasePaths.append(' ');
					}
					LOGGER.info(" * Reference genotype data: " + genotypeBasePaths);
					LOGGER.info(" * Reference genotype data type: " + genotypeType.getName());
				}
				LOGGER.info(" * Gene window extend in bases: " + windowExtend);
				LOGGER.info(" * Number of permutations: " + numberOfPermutations);
				LOGGER.info(" * Max correlation between variants: " + maxRBetweenVariants);
				LOGGER.info(" * Number of threads to use: " + numberOfThreadsToUse);
				LOGGER.info(" * Gene info file: " + geneInfoFile.getAbsolutePath());
				break;
		}
		
		LOGGER.info(" * Debug mode: " + (debugMode ? "on" : "off"));

	}

	public String[] getGenotypeBasePath() {
		return genotypeBasePath;
	}

	public RandomAccessGenotypeDataReaderFormats getGenotypeType() {
		return genotypeType;
	}

	public String getOutputBasePath() {
		return outputBasePath;
	}

	public File getLogFile() {
		return logFile;
	}

	public String getGwasZscoreMatrixPath() {
		return gwasZscoreMatrixPath;
	}

	public int getNumberOfPermutations() {
		return numberOfPermutations;
	}

	public int getWindowExtend() {
		return windowExtend;
	}

	public double getMaxRBetweenVariants() {
		return maxRBetweenVariants;
	}

	public static int getNumberOfThreadsToUse() {
		return numberOfThreadsToUse;
	}

	public File getGeneInfoFile() {
		return geneInfoFile;
	}

	public Depict2Mode getMode() {
		return mode;
	}

	public File getGenotypeSamplesFile() {
		return genotypeSamplesFile;
	}

}
