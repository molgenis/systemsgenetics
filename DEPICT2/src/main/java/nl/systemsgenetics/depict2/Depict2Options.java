/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import edu.emory.mathcs.utils.ConcurrencyUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;

import static nl.systemsgenetics.depict2.Depict2.LARGE_INT_FORMAT;

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
	private final int numberOfPermutations;
	private final long numberOfPermutationsRescue;
	private final int windowExtend;
	private final double maxRBetweenVariants;
	private final File logFile;
	private final boolean debugMode;
	private final File debugFolder;
	private final File intermediateFolder;
	private final boolean pvalueToZscore;
	private final List<PathwayDatabase> pathwayDatabases;
	private final File conversionColumnIncludeFilter;
	private final boolean correctForLambdaInflation;
	private final int permutationPathwayEnrichment;
	private final int permutationGeneCorrelations;
	private final boolean ignoreGeneCorrelations;
	private final double genePruningR;
	private final boolean forceNormalGenePvalues;
	private final boolean forceNormalPathwayPvalues;
	private final boolean normalizeEigenvectors;
	private final int geneCorrelationWindow;
	private final boolean excludeHla;
	private boolean corMatrixZscores = false;
	private String[] columnsToExtract = null; //Colums to extract when doing CONVERT_BIN or CONVERT_EXP
	private final File variantFilterFile;
	private boolean saveOuputAsExcelFiles;
	private final File variantGeneLinkingFile;
	private final boolean saveUsedVariantsPerGene;
	private final double mafFilter;

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
				+ "* CONVERT_TXT - Convert a txt z-score matrix to binary. Use --gwas, --output and optionally --pvalueToZscore if the matrix contains p-values instead of z-scores.\n"
				+ "* CONVERT_TXT_MERGE - Merge multiple txt pvalue files into one matrix containing only overlapping snps"
				+ "* CONVERT_BIN - Convert a binary matrix to a txt. Use --gwas and --output optionally --columnsToExtract\n"
				+ "* CONVERT_EXP - Convert a tab seperated expression matrix and normalize genes. Use --gwas (for exp data) and --output optionally --columnsToExtract\n"
				+ "* TRANSPOSE_BIN - Transposes a binary matrix. Use --gwas and --output\n"
				+ "* CONVERT_EQTL - Convert binary matrix with eQTL z-scores from our pipeline. Use --gwas and --output"
				+ "* CONVERT_GTEX - Convert Gtex median tissue GCT file. Use --gwas for the GCT file and --output"
				+ "* CORRELATE_GENES - Create gene correlation matrix with 0 on diagnonal. Use --gwas as input matrix (genes on row, tab sepperated), --output and --genes. Optionally use --corZscore to create Z-score matrix"
		);
		OptionBuilder.withLongOpt("mode");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("m"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("GWAS Z-sccore binary matrix. Rows variants, Cols phenotypes. Without .dat");
		OptionBuilder.withLongOpt("gwas");
		OPTIONS.addOption(OptionBuilder.create("g"));

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
				+ "* TRITYPER - TriTyper format folder");
		OptionBuilder.withLongOpt("referenceGenotypeFormat");
		OPTIONS.addOption(OptionBuilder.create("R"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("The output path");
		OptionBuilder.withLongOpt("output");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("o"));

		OptionBuilder.withArgName("int");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Maximum number of calculation threads");
		OptionBuilder.withLongOpt("threads");
		OPTIONS.addOption(OptionBuilder.create("t"));

		OptionBuilder.withArgName("int");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Number of initial permutations before using Farebrother's RUBEN algorithm to determine gene p-values. Recommended: 100,000; min: 10,000; max: " + LARGE_INT_FORMAT.format(GenePvalueCalculator.MAX_ROUND_1_RESCUE));
		OptionBuilder.withLongOpt("permutations");
		OPTIONS.addOption(OptionBuilder.create("p"));

		OptionBuilder.withArgName("int");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Number of permutations to use as fallback incase Farebrother's RUBEN algorithm failed. Optional but recommende to do atleast: 100,000,000. ");
		OptionBuilder.withLongOpt("permutationsRescue");
		OPTIONS.addOption(OptionBuilder.create("pr"));

		OptionBuilder.withArgName("int");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Number of initial permutations before using Farebrother's RUBEN algorithm to determine gene p-values ");
		OptionBuilder.withLongOpt("permutations");
		OPTIONS.addOption(OptionBuilder.create("p"));

		OptionBuilder.withArgName("int");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Number of bases to add left and right of gene window, use -1 to not have any variants in window");
		OptionBuilder.withLongOpt("window");
		OPTIONS.addOption(OptionBuilder.create("w"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Variant gene mapping. col 1: variant ID col 2: ENSG ID. (optional)");
		OptionBuilder.withLongOpt("variantGene");
		OPTIONS.addOption(OptionBuilder.create("vg"));

		OptionBuilder.withArgName("double");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Max correlation between variants to use (recommend = 0.95)");
		OptionBuilder.withLongOpt("variantCorrelation");
		OPTIONS.addOption(OptionBuilder.create("v"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
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
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Optional file with columns to select during conversion");
		OptionBuilder.withLongOpt("cols");
		OPTIONS.addOption(OptionBuilder.create("co"));

		OptionBuilder.withArgName("boolean");
		OptionBuilder.withDescription("Correct the GWAS for the lambda inflation");
		OptionBuilder.withLongOpt("correctLambda");
		OPTIONS.addOption(OptionBuilder.create("cl"));

		OptionBuilder.withArgName("boolean");
		OptionBuilder.withDescription("Ignore gene correlations in pathway enrichment (not recommended)");
		OptionBuilder.withLongOpt("ignoreGeneCorrelations");
		OPTIONS.addOption(OptionBuilder.create("igc"));

		OptionBuilder.withArgName("int");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Window in bases to calculate gene correlations for GLS of pathway enrichment");
		OptionBuilder.withLongOpt("geneCorrelationWindow");
		OPTIONS.addOption(OptionBuilder.create("gcw"));

		OptionBuilder.withArgName("int");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Number of random phenotypes to use to determine gene correlations");
		OptionBuilder.withLongOpt("permutationGeneCorrelations");
		OPTIONS.addOption(OptionBuilder.create("pgc"));

		OptionBuilder.withArgName("int");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Number of random phenotypes to use to determine null distribution pathway enrichment");
		OptionBuilder.withLongOpt("permutationPathwayEnrichment");
		OPTIONS.addOption(OptionBuilder.create("ppe"));

		OptionBuilder.withArgName("double");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Exclude correlated genes in pathway enrichments");
		OptionBuilder.withLongOpt("genePruningR");
		OPTIONS.addOption(OptionBuilder.create("gpr"));

		OptionBuilder.withArgName("boolean");
		OptionBuilder.withDescription("Force normal gene p-values before pathway enrichtment");
		OptionBuilder.withLongOpt("forceNormalGenePvalues");
		OPTIONS.addOption(OptionBuilder.create("fngp"));

		OptionBuilder.withArgName("boolean");
		OptionBuilder.withDescription("Force normal pathway p-values before pathway enrichtment");
		OptionBuilder.withLongOpt("forceNormalPathwayPvalues");
		OPTIONS.addOption(OptionBuilder.create("fnpp"));

		OptionBuilder.withArgName("boolean");
		OptionBuilder.withDescription("Exclude HLA region during pathway enrichments (chr6 20mb - 40mb)");
		OptionBuilder.withLongOpt("excludeHla");
		OPTIONS.addOption(OptionBuilder.create("eh"));

		OptionBuilder.withArgName("name=path");
		OptionBuilder.hasArgs();
		OptionBuilder.withValueSeparator();
		OptionBuilder.withDescription("Eigenvectors and princial componentens [path]_ev.txt & [path]_pc.txt");
		OptionBuilder.withLongOpt("pca");
		OPTIONS.addOption(OptionBuilder.create("pca"));

		OptionBuilder.withArgName("boolean");
		OptionBuilder.withDescription("When using --mode CORRELATE_GENES to correlate genes save results as Z-scores instead of r's.");
		OptionBuilder.withLongOpt("corZscore");
		OPTIONS.addOption(OptionBuilder.create("cz"));

		OptionBuilder.withArgName("boolean");
		OptionBuilder.withDescription("When using --mode CORRELATE_GENES first normalize the eigen vectors");
		OptionBuilder.withLongOpt("normalizeEigenvectors");
		OPTIONS.addOption(OptionBuilder.create("ne"));

		OptionBuilder.withArgName("boolean");
		OptionBuilder.withDescription("Save all the variants used per gene to calculate the gene p-value (Warning this will create a very large file)");
		OptionBuilder.withLongOpt("saveUsedVariantsPerGene");
		OPTIONS.addOption(OptionBuilder.create("uvg"));

		OptionBuilder.withArgName("strings");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Column names to extract when running --mode CONVERT_BIN or CONVERT_EXP");
		OptionBuilder.withLongOpt("columnsToExtract");
		OPTIONS.addOption(OptionBuilder.create("cte"));

		OptionBuilder.withArgName("strings");
		OptionBuilder.withDescription("Save enrichement results also as excel files. Will generate 1 file per input phenotype");
		OptionBuilder.withLongOpt("saveExcel");
		OPTIONS.addOption(OptionBuilder.create("se"));

		OptionBuilder.withArgName("double");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Minimum MAF");
		OptionBuilder.withLongOpt("maf");
		OPTIONS.addOption(OptionBuilder.create("maf"));

	}

	public Depict2Options(String... args) throws ParseException {

		final CommandLineParser parser = new PosixParser();
		final CommandLine commandLine = parser.parse(OPTIONS, args, false);

		if (commandLine.hasOption('t')) {
			try {
				numberOfThreadsToUse = Integer.parseInt(commandLine.getOptionValue('t'));
				System.setProperty("Djava.util.concurrent.ForkJoinPool.common.parallelism", commandLine.getOptionValue('t'));
				ConcurrencyUtils.setNumberOfThreads(numberOfThreadsToUse);
			} catch (NumberFormatException e) {
				throw new ParseException("Error parsing --threads \"" + commandLine.getOptionValue('t') + "\" is not an int");
			}
		}

		outputBasePath = new File(commandLine.getOptionValue('o'));
		logFile = new File(outputBasePath + ".log");
		debugMode = commandLine.hasOption('d');
		debugFolder = new File(outputBasePath + "_debugFiles");
		intermediateFolder = new File(outputBasePath + "_intermediates");
		ignoreGeneCorrelations = commandLine.hasOption("igc");
		correctForLambdaInflation = commandLine.hasOption("cl");
		forceNormalGenePvalues = commandLine.hasOption("fngp");
		forceNormalPathwayPvalues = commandLine.hasOption("fnpp");
		excludeHla = commandLine.hasOption("eh");
		normalizeEigenvectors = commandLine.hasOption("ne");
		saveOuputAsExcelFiles = commandLine.hasOption("se");
		saveUsedVariantsPerGene = commandLine.hasOption("uvg");

		try {
			mode = Depict2Mode.valueOf(commandLine.getOptionValue("m").toUpperCase());
		} catch (IllegalArgumentException e) {
			throw new ParseException("Error parsing --mode \"" + commandLine.getOptionValue("m") + "\" is not a valid mode");
		}

		if (mode == Depict2Mode.CONVERT_TXT || mode == Depict2Mode.CONVERT_TXT_MERGE || mode == Depict2Mode.RUN || mode == Depict2Mode.CONVERT_EQTL || mode == Depict2Mode.FIRST1000 || mode == Depict2Mode.CONVERT_GTEX || mode == Depict2Mode.CONVERT_BIN || mode == Depict2Mode.SPECIAL || mode == Depict2Mode.CORRELATE_GENES || mode == Depict2Mode.TRANSPOSE || mode == Depict2Mode.CONVERT_EXP) {

			if (!commandLine.hasOption("g")) {
				throw new ParseException("Please provide --gwas for mode: " + mode.name());
			}

			gwasZscoreMatrixPath = new File(commandLine.getOptionValue('g'));
		} else {
			gwasZscoreMatrixPath = null;
		}

		if (mode == Depict2Mode.CONVERT_TXT || mode == Depict2Mode.CONVERT_TXT_MERGE) {
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

		switch (mode) {
			case RUN:
			case RUN2:
				if (!commandLine.hasOption("pgc")) {
					throw new ParseException("--permutationGeneCorrelations not specified");
				} else {
					try {
						permutationGeneCorrelations = Integer.parseInt(commandLine.getOptionValue("pgc"));
					} catch (NumberFormatException e) {
						throw new ParseException("Error parsing --permutationGeneCorrelations \"" + commandLine.getOptionValue("pgc") + "\" is not an int");
					}
				}
				if (!commandLine.hasOption("ppe")) {
					throw new ParseException("--permutationPathwayEnrichment not specified");
				} else {
					try {
						permutationPathwayEnrichment = Integer.parseInt(commandLine.getOptionValue("ppe"));
					} catch (NumberFormatException e) {
						throw new ParseException("Error parsing --permutationPathwayEnrichment \"" + commandLine.getOptionValue("ppe") + "\" is not an int");
					}
				}
				if (!commandLine.hasOption("gpr")) {
					throw new ParseException("--genePruningR not specified");
				} else {
					try {
						genePruningR = Double.parseDouble(commandLine.getOptionValue("gpr"));
					} catch (NumberFormatException e) {
						throw new ParseException("Error parsing --genePruningR \"" + commandLine.getOptionValue("gpr") + "\" is not an double");
					}
				}
				if (!commandLine.hasOption("ge")) {
					throw new ParseException("--genes not specified");
				} else {
					geneInfoFile = new File(commandLine.getOptionValue("ge"));
				}
				if (!commandLine.hasOption("gcw")) {
					throw new ParseException("--geneCorrelationWindow not specified");
				} else {
					try {
						geneCorrelationWindow = Integer.parseInt(commandLine.getOptionValue("gcw"));
					} catch (NumberFormatException e) {
						throw new ParseException("Error parsing --geneCorrelationWindow \"" + commandLine.getOptionValue("gcw") + "\" is not an int");
					}
				}
				pathwayDatabases = parsePd(commandLine);
				break;
			case CORRELATE_GENES:
				if (commandLine.hasOption("ge")) {
					geneInfoFile = new File(commandLine.getOptionValue("ge"));
				} else {
					geneInfoFile = null;
				}
				corMatrixZscores = commandLine.hasOption("cz");
				pathwayDatabases = null;
				permutationGeneCorrelations = 0;
				permutationPathwayEnrichment = 0;
				genePruningR = 0;
				geneCorrelationWindow = 0;
				break;
			case CONVERT_BIN:
			case CONVERT_EXP:
				if (commandLine.hasOption("cte")) {
					columnsToExtract = commandLine.getOptionValues("cte");
				}
			default:
				pathwayDatabases = null;
				permutationGeneCorrelations = 0;
				permutationPathwayEnrichment = 0;
				genePruningR = 0;
				geneInfoFile = null;
				geneCorrelationWindow = 0;
				break;
		}

		switch (mode) {
			case RUN:

				//variantFilter
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

				if (!commandLine.hasOption("p")) {
					throw new ParseException("--permutations not specified");
				} else {

					try {
						numberOfPermutations = Integer.parseInt(commandLine.getOptionValue('p'));
					} catch (NumberFormatException e) {
						throw new ParseException("Error parsing --permutations \"" + commandLine.getOptionValue('p') + "\" is not an int");
					}

					if (numberOfPermutations > GenePvalueCalculator.MAX_ROUND_1_RESCUE) {
						throw new ParseException("Error parsing --permutations max is: " + GenePvalueCalculator.MAX_ROUND_1_RESCUE);
					}

				}
				if (!commandLine.hasOption("pr")) {
					numberOfPermutationsRescue = numberOfPermutations;
				} else {

					try {
						numberOfPermutationsRescue = Long.parseLong(commandLine.getOptionValue("pr"));
					} catch (NumberFormatException e) {
						throw new ParseException("Error parsing --permutationsRescue \"" + commandLine.getOptionValue('p') + "\" is not an long");
					}

					if (numberOfPermutationsRescue < numberOfPermutations) {
						throw new ParseException("--permutationsRescue must be atleast equal to --permutations. If no extra permutations are wanted this option can be omited");
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

				if (commandLine.hasOption("vg")) {
					variantGeneLinkingFile = new File(commandLine.getOptionValue("vg"));
				} else {
					variantGeneLinkingFile = null;
				}

				if (windowExtend < 0 && variantGeneLinkingFile == null) {
					throw new ParseException("No window specified but also no variant gene linking.");
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
				break;
			case RUN2:
				if (pathwayDatabases.isEmpty()) {
					throw new ParseException("The option --pathwayDatabase is needed for mode=RUN2");
				}
				genotypeBasePath = null;
				genotypeType = null;
				genotypeSamplesFile = null;
				maxRBetweenVariants = 0d;
				numberOfPermutations = 0;
				numberOfPermutationsRescue = 0;
				windowExtend = 0;
				variantFilterFile = null;
				variantGeneLinkingFile = null;
				mafFilter = 0;
				break;
			default:
				genotypeBasePath = null;
				genotypeType = null;
				genotypeSamplesFile = null;
				maxRBetweenVariants = 0d;
				numberOfPermutations = 0;
				numberOfPermutationsRescue = 0;
				windowExtend = 0;
				variantFilterFile = null;
				variantGeneLinkingFile = null;
				mafFilter = 0;
				break;
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
				if (conversionColumnIncludeFilter != null) {
					LOGGER.info(" * Columns to include: " + conversionColumnIncludeFilter.getAbsolutePath());
				}
				LOGGER.info(" * Convert p-values to Z-score: " + (pvalueToZscore ? "on" : "off"));
				break;
			case CONVERT_TXT_MERGE:
				LOGGER.info(" * Gwas Z-score matrix paths: " + gwasZscoreMatrixPath.getAbsolutePath());
				LOGGER.info(" * Convert p-values to Z-score: " + (pvalueToZscore ? "on" : "off"));
				break;
			case CONVERT_BIN:
				LOGGER.info(" * Gwas Z-score matrix: " + gwasZscoreMatrixPath.getAbsolutePath());
				if (columnsToExtract != null) {
					LOGGER.info(" * Columns to extract: " + String.join(" ", columnsToExtract));
				}
				break;
			case CONVERT_EXP:
				LOGGER.info(" * Expression matrix: " + gwasZscoreMatrixPath.getAbsolutePath());
				if (columnsToExtract != null) {
					LOGGER.info(" * Columns to extract: " + String.join(" ", columnsToExtract));
				}
				break;
			case TRANSPOSE:
				LOGGER.info(" * Matrix to transpose: " + gwasZscoreMatrixPath.getAbsolutePath());
				break;
			case CORRELATE_GENES:
				LOGGER.info(" * Gwas Z-score matrix: " + gwasZscoreMatrixPath.getAbsolutePath());
				LOGGER.info(" * First normalize eigen vectors: " + (normalizeEigenvectors ? "on" : "off"));
				LOGGER.info(" * Convert r to Z-score: " + (corMatrixZscores ? "on" : "off"));
				if (geneInfoFile != null) {
					LOGGER.info(" * Genes to include file: " + geneInfoFile.getAbsolutePath());
				}
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
				if (mafFilter != 0) {
					LOGGER.info(" * MAF filter: " + mafFilter);
				}
				LOGGER.info(" * Gene window extend in bases: " + (windowExtend < 0 ? "Not selecting variants in a window" : LARGE_INT_FORMAT.format(windowExtend)));
				if (variantGeneLinkingFile != null) {
					LOGGER.info(" * Variant to gene linking file: " + variantGeneLinkingFile);
				}
				LOGGER.info(" * Initial number of permutations to calculate gene p-values: " + LARGE_INT_FORMAT.format(numberOfPermutations));
				LOGGER.info(" * Max number of rescue permutations to calculate gene p-values if RUBEN has failed: " + LARGE_INT_FORMAT.format(numberOfPermutationsRescue));
				LOGGER.info(" * Max correlation between variants: " + maxRBetweenVariants);
				LOGGER.info(" * Correcting for lambda inflation: " + (correctForLambdaInflation ? "on" : "off"));
				LOGGER.info(" * Save which variants that are used per gene to calculate the gene p-value: " + (saveUsedVariantsPerGene ? "on" : "off"));
				if (variantFilterFile != null) {
					LOGGER.info(" * Confining analysis to variants in this file: " + variantFilterFile.getAbsolutePath());
				}
				logSharedRun1Run2();

				break;
			case RUN2:
				logSharedRun1Run2();

				break;

		}

		LOGGER.info(" * Debug mode: " + (debugMode ? "on (this will result in many intermediate output files)" : "off"));

	}

	private void logSharedRun1Run2() {

		if (pvalueToZscore) {
			LOGGER.info("WARNING --pvalueToZscore is set but only effective for mode: CONVERT_TXT");
		}
		LOGGER.info(" * Gene info file: " + geneInfoFile.getAbsolutePath());
		LOGGER.info(" * Number of threads to use: " + numberOfThreadsToUse);

		LOGGER.info(" * Number of permutations to use to calculate gene correlations: " + LARGE_INT_FORMAT.format(permutationGeneCorrelations));
		LOGGER.info(" * Number of permutations to use for pathway enrichments: " + LARGE_INT_FORMAT.format(permutationPathwayEnrichment));
		LOGGER.info(" * Window to calculate gene correlations for GLS: " + LARGE_INT_FORMAT.format(geneCorrelationWindow));
		LOGGER.info(" * Gene pruning r: " + genePruningR);
		LOGGER.info(" * Ignoring gene correlations: " + (ignoreGeneCorrelations ? "on" : "off"));
		LOGGER.info(" * Force normal gene p-values: " + (forceNormalGenePvalues ? "on" : "off"));
		LOGGER.info(" * Force normal pathway p-values: " + (forceNormalPathwayPvalues ? "on" : "off"));
		LOGGER.info(" * Exclude HLA during enrichment analysis: " + (excludeHla ? "on" : "off"));
		LOGGER.info(" * Save output as excel files: " + (saveOuputAsExcelFiles ? "on" : "off"));
		logPathwayDatabases();

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

	public int getPermutationPathwayEnrichment() {
		return permutationPathwayEnrichment;
	}

	public int getPermutationGeneCorrelations() {
		return permutationGeneCorrelations;
	}

	public boolean isIgnoreGeneCorrelations() {
		return ignoreGeneCorrelations;
	}

	public double getGenePruningR() {
		return genePruningR;
	}

	public boolean isForceNormalGenePvalues() {
		return forceNormalGenePvalues;
	}

	public boolean isForceNormalPathwayPvalues() {
		return forceNormalPathwayPvalues;
	}

	public int getGeneCorrelationWindow() {
		return geneCorrelationWindow;
	}

	public boolean isExcludeHla() {
		return excludeHla;
	}

	public long getNumberOfPermutationsRescue() {
		return numberOfPermutationsRescue;
	}

	public boolean isCorMatrixZscores() {
		return corMatrixZscores;
	}

	public String[] getColumnsToExtract() {
		return columnsToExtract;
	}

	public File getDebugFolder() {
		return debugFolder;
	}

	public File getIntermediateFolder() {
		return intermediateFolder;
	}

	public boolean isNormalizeEigenvectors() {
		return normalizeEigenvectors;
	}

	public File getVariantFilterFile() {
		return variantFilterFile;
	}

	public boolean isSaveOuputAsExcelFiles() {
		return saveOuputAsExcelFiles;
	}

	public File getVariantGeneLinkingFile() {
		return variantGeneLinkingFile;
	}

	public boolean isSaveUsedVariantsPerGene() {
		return saveUsedVariantsPerGene;
	}

	public double getMafFilter() {
		return mafFilter;
	}

}
