/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
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
	private static int numberOfThreadsToUse = Runtime.getRuntime().availableProcessors();//Might be changed
	private static final Logger LOGGER = Logger.getLogger(Depict2Options.class);

	private final Depict2Mode mode;

	private final String[] genotypeBasePath;
	private final RandomAccessGenotypeDataReaderFormats genotypeType;
	private final File genotypeSamplesFile;
	private final File outputBasePath;
	private final File geneInfoFile;
	private final File gwasZscoreMatrixPath;
	private final long numberOfPermutations;
	private final int windowExtend;
	private final double maxRBetweenVariants;
	private final File logFile;
	private final boolean debugMode;
	private final boolean pvalueToZscore;
	private final List<PathwayDatabase> pathwayDatabases;
	private final File conversionColumnIncludeFilter;
	private final boolean correctForLambdaInflation;

	public boolean isDebugMode() {
		return debugMode;
	}

	static {

		OPTIONS = new Options();

		OptionBuilder.withArgName("string");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("On of the following modes:\n"
				+ "* RUN - Run the DEPICT2 prioritization.\n"
				+ "* RUN2 - Run the DEPICT2 prioritization starting at stage 2.\n"
				+ "* RUN3 - Run the DEPICT2 prioritization starting at stage 3.\n"
				+ "* CONVERT_TXT - Convert a txt z-score matrix to binary. Use --gwas, --output and optionally --pvalueToZscore if the matrix contains p-values instead of z-scores.\n"
				+ "* CONVERT_TXT - Convert a binary matrix to a txt. Use --gwas and --output\n"
				+ "* CONVERT_EQTL - Convert binary matrix with eQTL z-scores from our pipeline. Use --gwas and --output"
				+ "* CONVERT_GTEX - Convert Gtex median tissue GCT file. Use --gwas for the GCT file and --output");
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
		OptionBuilder.isRequired();
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
		OptionBuilder.withDescription("File with gene info. col1: geneName (ensg) col2: chr col3: startPos col4: stopPos col5: geneType col6: chrArm");
		OptionBuilder.withLongOpt("genes");
		OPTIONS.addOption(OptionBuilder.create("ge"));

		OptionBuilder.withArgName("boolean");
		OptionBuilder.withDescription("Activate debug mode. This will result in a more verbose log file and will save many intermediate results to files. Not recommended for large analysis.");
		OptionBuilder.withLongOpt("debug");
		OPTIONS.addOption(OptionBuilder.create("d"));

		OptionBuilder.withArgName("boolean");
		OptionBuilder.withDescription("When mode=CONVERT_TXT convert p-values to z-scores");
		OptionBuilder.withLongOpt("pvalueToZscore");
		OPTIONS.addOption(OptionBuilder.create("p2z"));

		OptionBuilder.withArgName("name=path");
		OptionBuilder.hasArgs();
		OptionBuilder.withValueSeparator();
		OptionBuilder.withDescription("Pathway databases, binary matrix with p-values for genes");
		OptionBuilder.withLongOpt("pathwayDatabase");
		OPTIONS.addOption(OptionBuilder.create("pd"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Optional file with columns to select during conversion");
		OptionBuilder.withLongOpt("cols");
		OPTIONS.addOption(OptionBuilder.create("co"));

		OptionBuilder.withArgName("boolean");
		OptionBuilder.withDescription("Correct the GWAS for the lambda inflation");
		OptionBuilder.withLongOpt("correctLambda");
		OPTIONS.addOption(OptionBuilder.create("cl"));

	}

	public Depict2Options(String... args) throws ParseException {

		final CommandLineParser parser = new PosixParser();
		final CommandLine commandLine = parser.parse(OPTIONS, args, false);

		if (commandLine.hasOption('t')) {
			try {
				numberOfThreadsToUse = Integer.parseInt(commandLine.getOptionValue('t'));
				System.setProperty("Djava.util.concurrent.ForkJoinPool.common.parallelism", commandLine.getOptionValue('t'));
			} catch (NumberFormatException e) {
				throw new ParseException("Error parsing --threads \"" + commandLine.getOptionValue('t') + "\" is not an int");
			}
		}

		outputBasePath = new File(commandLine.getOptionValue('o'));
		logFile = new File(outputBasePath + ".log");
		debugMode = commandLine.hasOption('d');
		correctForLambdaInflation = commandLine.hasOption("cl");

		try {
			mode = Depict2Mode.valueOf(commandLine.getOptionValue("m").toUpperCase());
		} catch (IllegalArgumentException e) {
			throw new ParseException("Error parsing --mode \"" + commandLine.getOptionValue("m") + "\" is not a valid mode");
		}

		if (mode == Depict2Mode.CONVERT_TXT || mode == Depict2Mode.RUN || mode == Depict2Mode.CONVERT_EQTL || mode == Depict2Mode.FIRST1000 || mode == Depict2Mode.CONVERT_GTEX || mode == Depict2Mode.CONVERT_BIN || mode == Depict2Mode.SPECIAL) {

			if (!commandLine.hasOption("g")) {
				throw new ParseException("Please provide --gwas for mode: " + mode.name());
			}

			gwasZscoreMatrixPath = new File(commandLine.getOptionValue('g'));
		} else {
			gwasZscoreMatrixPath = null;
		}

		if (mode == Depict2Mode.CONVERT_TXT) {
			pvalueToZscore = commandLine.hasOption("p2z");
			if (commandLine.hasOption("co")) {
				conversionColumnIncludeFilter = new File(commandLine.getOptionValue("co"));
			} else {
				conversionColumnIncludeFilter = null;
			}
		} else {
			pvalueToZscore = false;
			conversionColumnIncludeFilter = null;

		}

		if (mode == Depict2Mode.RUN) {

			if (!commandLine.hasOption("p")) {
				throw new ParseException("--permutations not specified");
			} else {

				try {
					numberOfPermutations = Long.parseLong(commandLine.getOptionValue('p'));
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

			pathwayDatabases = parsePd(commandLine);

		} else if (mode == Depict2Mode.RUN2 || mode == Depict2Mode.RUN3) {

			if (!commandLine.hasOption("ge")) {
				throw new ParseException("--genes not specified");
			} else {
				geneInfoFile = new File(commandLine.getOptionValue("ge"));
			}

			pathwayDatabases = parsePd(commandLine);

			if (mode == Depict2Mode.RUN3 && pathwayDatabases.isEmpty()) {
				throw new ParseException("The option --pathwayDatabase is needed for mode=RUN3");
			}

			genotypeBasePath = null;
			genotypeType = null;
			genotypeSamplesFile = null;
			maxRBetweenVariants = 0d;
			numberOfPermutations = 0;
			windowExtend = 0;

		} else {
			pathwayDatabases = null;
			genotypeBasePath = null;
			genotypeType = null;
			genotypeSamplesFile = null;
			geneInfoFile = null;
			maxRBetweenVariants = 0d;
			numberOfPermutations = 0;
			windowExtend = 0;
		}

	}

	private List<PathwayDatabase> parsePd(final CommandLine commandLine) throws ParseException {

		final List<PathwayDatabase> pathwayDatabasesTmp;

		if (commandLine.hasOption("pd")) {

			String[] pdValues = commandLine.getOptionValues("pd");

			if (pdValues.length % 2 != 0) {
				throw new ParseException("Error parsing --pathwayDatabase. Must be in name=database format");
			}

			final HashSet<String> duplicateChecker = new HashSet<>();
			pathwayDatabasesTmp = new ArrayList<>();

			for (int i = 0; i < pdValues.length; i += 2) {

				if (!duplicateChecker.add(pdValues[i])) {
					throw new ParseException("Error parsing --pathwayDatabase. Duplicate database name found");
				}

				pathwayDatabasesTmp.add(new PathwayDatabase(pdValues[i], pdValues[i + 1]));

			}

		} else {

			pathwayDatabasesTmp = Collections.emptyList();

		}

		return pathwayDatabasesTmp;
	}

	public static void printHelp() {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp(" ", OPTIONS);
	}

	public void printOptions() {

		LOGGER.info("Supplied options:");

		LOGGER.info(" * Mode: " + mode.name());
		LOGGER.info(" * Ouput path: " + outputBasePath.getAbsolutePath());

		switch (mode) {
			case CONVERT_EQTL:
				LOGGER.info(" * eQTL Z-score matrix: " + gwasZscoreMatrixPath.getAbsolutePath());
				if (pvalueToZscore) {
					LOGGER.info("WARNING --pvalueToZscore is set but only effective for mode: CONVERT_TXT");
				}
				break;
			case CONVERT_GTEX:
				LOGGER.info(" * Gtex median tissue expression GCT file: " + gwasZscoreMatrixPath.getAbsolutePath());
				if (pvalueToZscore) {
					LOGGER.info("WARNING --pvalueToZscore is set but only effective for mode: CONVERT_TXT");
				}
				break;
			case CONVERT_TXT:
				LOGGER.info(" * Gwas Z-score matrix: " + gwasZscoreMatrixPath.getAbsolutePath());
				LOGGER.info(" * Columns to include: " + conversionColumnIncludeFilter.getAbsolutePath());
				LOGGER.info(" * Convert p-values to Z-score: " + (pvalueToZscore ? "on" : "off"));
				break;
			case CONVERT_BIN:
				LOGGER.info(" * Gwas Z-score matrix: " + gwasZscoreMatrixPath.getAbsolutePath());
				break;
			case RUN:
				LOGGER.info(" * Gwas Z-score matrix: " + gwasZscoreMatrixPath.getAbsolutePath());

				if (genotypeBasePath != null) {
					StringBuilder genotypeBasePaths = new StringBuilder();
					for (String path : genotypeBasePath) {
						genotypeBasePaths.append(new File(path).getAbsolutePath());
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
				if (pvalueToZscore) {
					LOGGER.info("WARNING --pvalueToZscore is set but only effective for mode: CONVERT_TXT");
				}
				LOGGER.info(" * Correcting for lambda inflation: " + (correctForLambdaInflation ? "on" : "off"));

				logPathwayDatabases();

				break;
			case RUN2:
			case RUN3:
				LOGGER.info(" * Gene info file: " + geneInfoFile.getAbsolutePath());
				logPathwayDatabases();
				break;

		}

		LOGGER.info(" * Debug mode: " + (debugMode ? "on (this will result in many intermediate output files)" : "off"));

	}

	private void logPathwayDatabases() {
		if (pathwayDatabases.size() > 0) {
			LOGGER.info(" * The following pathway databases have been specified:");
			for (PathwayDatabase database : pathwayDatabases) {
				LOGGER.info("    - " + database.getName() + " at: " + database.getLocation());
			}
		} else {
			LOGGER.info(" * No pathway databases specified");
		}
	}

	public String[] getGenotypeBasePath() {
		return genotypeBasePath;
	}

	public RandomAccessGenotypeDataReaderFormats getGenotypeType() {
		return genotypeType;
	}

	public String getOutputBasePath() {
		return outputBasePath.getPath();
	}

	public File getLogFile() {
		return logFile;
	}

	public String getGwasZscoreMatrixPath() {
		return gwasZscoreMatrixPath.getPath();
	}

	public long getNumberOfPermutations() {
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

	public boolean isPvalueToZscore() {
		return pvalueToZscore;
	}

	public List<PathwayDatabase> getPathwayDatabases() {
		return pathwayDatabases;
	}

	public File getConversionColumnIncludeFilter() {
		return conversionColumnIncludeFilter;
	}

	public boolean correctForLambdaInflation() {
		return correctForLambdaInflation;
	}

}
