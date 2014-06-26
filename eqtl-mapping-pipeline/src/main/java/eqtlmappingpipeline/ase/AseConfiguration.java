/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.ase;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import org.apache.commons.cli.CommandLine;
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
 * @author Patrick Deelen
 */
public class AseConfiguration {

	private static final Logger LOGGER;
	private static final Options OPTIONS;
	public static final String ENCODING = "UTF-8";
	private final List<File> inputFiles;
	private final File outputFolder;
	private final int minTotalReads;
	private final int minAlleleReads;
	private final double minAlleleReadFraction;
	private final File logFile;
	private final boolean debugMode;
	private final int minSamples;
	private final int threads;
	private final String[] refBasePaths;
	private final RandomAccessGenotypeDataReaderFormats refDataType;
	private final int refDataCacheSize;
	private final int maxTotalReads;
	private final File gtf;
	private final File sampleToRefSampleFile;

	static {

		LOGGER = Logger.getLogger(AseConfiguration.class);

		OPTIONS = new Options();

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Paths to one or more vcf.gz files or folders with each per chr 1 vcf.gz file. Tabix file must be present");
		OptionBuilder.withLongOpt("input");
		OPTIONS.addOption(OptionBuilder.create('i'));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("One or more text files with on each line a path to a vcf.gz file or a folder with per chr 1 vcf.gz file. Tabix file must be present. Can be combined with --input");
		OptionBuilder.withLongOpt("inputList");
		OPTIONS.addOption(OptionBuilder.create('l'));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Path to output folder");
		OptionBuilder.withLongOpt("output");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create('o'));

		OptionBuilder.withArgName("int");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Min number total reads for genotype");
		OptionBuilder.withLongOpt("minReads");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create('r'));

		OptionBuilder.withArgName("int");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Min number of reads per allele");
		OptionBuilder.withLongOpt("minAlleleReads");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create('a'));

		OptionBuilder.withArgName("int");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Min percentage of read per allele. Default 0. Can be used in combination --minAlleleReads");
		OptionBuilder.withLongOpt("minAllelePercentage");
		OPTIONS.addOption(OptionBuilder.create('p'));

		OptionBuilder.withArgName("int");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Min number of samples per ASE effect");
		OptionBuilder.withLongOpt("minNumSamples");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create('s'));

		OptionBuilder.withArgName("int");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Maximum number of threads to start for parallel reading of multiple VCF files. Defaults to number of cores");
		OptionBuilder.withLongOpt("threads");
		OPTIONS.addOption(OptionBuilder.create('t'));

		OptionBuilder.withArgName("basePath");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("The path to the reference genotypes. These genotypes will be used to determine if a sample is hetrozygous");
		OptionBuilder.withLongOpt("genotypes");
		OPTIONS.addOption(OptionBuilder.create("g"));

		OptionBuilder.withArgName("type");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("The input data type. If not defined will attempt to automatically select the first matching dataset on the specified path\n"
				+ "* PED_MAP - plink PED MAP files.\n"
				+ "* PLINK_BED - plink BED BIM FAM files.\n"
				+ "* VCF - bgziped vcf with tabix index file\n"
				+ "* VCFFOLDER - matches all bgziped vcf files + tabix index in a folder\n"
				+ "* SHAPEIT2 - shapeit2 phased haplotypes .haps & .sample\n"
				+ "* GEN - Oxford .gen & .sample\n"
				+ "* TRITYPER - TriTyper format folder");
		OptionBuilder.withLongOpt("genotypesType");
		OPTIONS.addOption(OptionBuilder.create("G"));

		OptionBuilder.withArgName("int");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Reference genotype data cache. Trade memory usage for speed");
		OptionBuilder.withLongOpt("cache");
		OPTIONS.addOption(OptionBuilder.create("c"));
		
		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Tab separated file with column 1 sample ID and column 2 sample ID in reference genotype data, no header. If only contains a mapping for a subset of samples the original identifier in the reference is used");
		OptionBuilder.withLongOpt("sampleCoupling");
		OPTIONS.addOption(OptionBuilder.create("sc"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Optional .gtf file for annotations of ASE effects. Must be grouped by chromosome");
		OptionBuilder.withLongOpt("gtf");
		OPTIONS.addOption(OptionBuilder.create("f"));

		OptionBuilder.withArgName("int");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Maximum number of total reads");
		OptionBuilder.withLongOpt("maxReads");
		OPTIONS.addOption(OptionBuilder.create("m"));

		OptionBuilder.withArgName("boolean");
		OptionBuilder.withDescription("Activate debug mode. This will result in a more verbose log file");
		OptionBuilder.withLongOpt("debug");
		OPTIONS.addOption(OptionBuilder.create('d'));

	}

	public AseConfiguration(String... args) throws ParseException {

		final CommandLine commandLine = new PosixParser().parse(OPTIONS, args, true);

		if (!commandLine.hasOption('i') && !commandLine.hasOption('l')) {
			throw new ParseException("At least --input or --inputList need to be supplied");
		}

		ArrayList<File> inputFilesTmp = new ArrayList<File>();

		if (commandLine.hasOption('i')) {
			String[] inputPaths = commandLine.getOptionValues('i');
			for (String inputPath : inputPaths) {
				inputFilesTmp.add(new File(inputPath));
			}
		}

		if (commandLine.hasOption('l')) {
			String[] inputListFilePaths = commandLine.getOptionValues('l');
			for (String inputListFilePath : inputListFilePaths) {

				try {
					BufferedReader inputListFileReader = new BufferedReader(new InputStreamReader(new FileInputStream(inputListFilePath), ENCODING));

					String line;
					while ((line = inputListFileReader.readLine()) != null) {
						inputFilesTmp.add(new File(line));
					}


				} catch (FileNotFoundException ex) {
					throw new ParseException("Could not find file with list of input files: " + inputListFilePath);
				} catch (UnsupportedEncodingException ex) {
					throw new RuntimeException("Fatal error", ex);
				} catch (IOException ex) {
					throw new ParseException("Error reading list of input files: " + inputListFilePath);
				}


			}
		}


		inputFiles = Collections.unmodifiableList(inputFilesTmp);

		outputFolder = new File(commandLine.getOptionValue('o'));

		logFile = new File(outputFolder, "ase.log");

		try {
			minTotalReads = Integer.parseInt(commandLine.getOptionValue('r'));
			if (minTotalReads <= 0) {
				throw new ParseException("--minReads must be larger than 0");
			}
		} catch (NumberFormatException e) {
			throw new ParseException("Error parsing --minReads \"" + commandLine.getOptionValue('r') + "\" is not an int");
		}

		try {
			minAlleleReads = Integer.parseInt(commandLine.getOptionValue('a'));
			if (minAlleleReads < 0) {
				throw new ParseException("--minAlleleReads must be positive");
			}
		} catch (NumberFormatException e) {
			throw new ParseException("Error parsing --minAlleleReads \"" + commandLine.getOptionValue('a') + "\" is not an int");
		}

		if (commandLine.hasOption('p')) {
			try {
				minAlleleReadFraction = Double.parseDouble(commandLine.getOptionValue('p')) / 100;
				if(minAlleleReadFraction < 0 || minAlleleReadFraction >= 0.5){
					throw new ParseException("--minAllelePercentage must be in interval [0, 50)");
				}
			} catch (NumberFormatException e) {
				throw new ParseException("Error parsing --minAllelePercentage \"" + commandLine.getOptionValue('p') + "\" is not a double");
			}
		} else {
			minAlleleReadFraction = 0;
		}

		try {
			minSamples = Integer.parseInt(commandLine.getOptionValue('s'));
			if (minSamples <= 0) {
				throw new ParseException("--minNumSamples must be larger than 0");
			}
		} catch (NumberFormatException e) {
			throw new ParseException("Error parsing --minNumSamples \"" + commandLine.getOptionValue('s') + "\" is not an int");
		}

		int availCores = Runtime.getRuntime().availableProcessors();
		if (commandLine.hasOption('t')) {
			try {
				int threadOption = Integer.parseInt(commandLine.getOptionValue('t'));
				if (threadOption <= 0) {
					throw new ParseException("--threads must be larger than 0");
				}
				threads = threadOption > availCores ? availCores : threadOption;
			} catch (NumberFormatException e) {
				throw new ParseException("Error parsing --threads \"" + commandLine.getOptionValue('t') + "\" is not an int");
			}
		} else {
			threads = availCores;
		}


		if (commandLine.hasOption('c')) {
			try {
				refDataCacheSize = Integer.parseInt(commandLine.getOptionValue('c'));
				if (refDataCacheSize < 0) {
					throw new ParseException("--cache must be positive");
				}
			} catch (NumberFormatException e) {
				throw new ParseException("Error parsing --cache \"" + commandLine.getOptionValue('c') + "\" is not an int");
			}
		} else {
			refDataCacheSize = 1000;
		}

		if (commandLine.hasOption('g')) {
			refBasePaths = commandLine.getOptionValues('g');
			try {
				if (commandLine.hasOption('G')) {
					refDataType = RandomAccessGenotypeDataReaderFormats.valueOf(commandLine.getOptionValue('G').toUpperCase());
				} else {
					if (refBasePaths[0].endsWith(".vcf")) {
						throw new ParseException("Only vcf.gz is supported. Please see manual on how to do create a vcf.gz file.");
					}
					try {
						refDataType = RandomAccessGenotypeDataReaderFormats.matchFormatToPath(refBasePaths[0]);
					} catch (GenotypeDataException e) {
						throw new ParseException("Unable to determine reference data type based on specified path. Please specify --genotypesType");
					}
				}
			} catch (IllegalArgumentException e) {
				throw new ParseException("Error parsing --genotypesType \"" + commandLine.getOptionValue('G') + "\" is not a valid input data format");
			}
		} else {
			refBasePaths = null;
			refDataType = null;
		}

		if (commandLine.hasOption('m')) {
			try {
				maxTotalReads = Integer.parseInt(commandLine.getOptionValue('m'));
				if (maxTotalReads <= 0) {
					throw new ParseException("--maxReads must be larger than 0");
				}
			} catch (NumberFormatException e) {
				throw new ParseException("Error parsing --maxReads \"" + commandLine.getOptionValue('m') + "\" is not an int");
			}
		} else {
			maxTotalReads = Integer.MAX_VALUE;
		}

		if (commandLine.hasOption('f')) {
			String gencodeGtfPath = commandLine.getOptionValue('f');
			gtf = new File(gencodeGtfPath);
		} else {
			gtf = null;
		}
		
		
		if (commandLine.hasOption("sc")) {
			sampleToRefSampleFile = new File(commandLine.getOptionValue("sc"));
		} else {
			sampleToRefSampleFile = null;
		}

		debugMode = commandLine.hasOption('d');

	}

	public void printOptions() {

		System.out.println("Interpreted arguments: ");


		System.out.println(" - Input files or folders (" + Ase.DEFAULT_NUMBER_FORMATTER.format(inputFiles.size()) + " in total): ");
		if (inputFiles.size() > 5) {
			System.out.println("  * Reading more than 5 files. See log file for list.");
		} else {
			for (File inputFile : inputFiles) {
				System.out.println("  * " + inputFile.getAbsolutePath());
			}
		}

		LOGGER.info("Input files or folders (" + Ase.DEFAULT_NUMBER_FORMATTER.format(inputFiles.size()) + " in total): ");
		for (File inputFile : inputFiles) {
			LOGGER.info(" * " + inputFile.getAbsolutePath());
		}

		System.out.println(" - Output folder: " + outputFolder.getAbsolutePath());
		LOGGER.info("Output folder: " + outputFolder.getAbsolutePath());

		System.out.println(" - Minimum number of reads per genotype: " + minTotalReads);
		LOGGER.info("Minimum number of reads per genotype: " + minTotalReads);

		System.out.println(" - Minimum number of reads per allele: " + minAlleleReads);
		LOGGER.info("Minimum number of reads per allele: " + minAlleleReads);

		System.out.println(" - Minimum percentage of reads per allele: " + minAlleleReadFraction * 100 + "%");
		LOGGER.info("Minimum percentage of reads per allele: " + minAlleleReadFraction * 100 + "%");

		if (maxTotalReads != Integer.MAX_VALUE) {
			System.out.println(" - Maximum number of reads per genotype: " + maxTotalReads);
			LOGGER.info("Maximum number of reads per genotype: " + maxTotalReads);
		}

		System.out.println(" - Minimum number of samples per ASE effect: " + minSamples);
		LOGGER.info("Minimum number of samples per ASE effect: " + minSamples);

		System.out.println(" - Number of threads to use: " + threads);
		LOGGER.info("Number of threads to use: " + threads);

		if (isRefSet()) {
			System.out.print(" - Reference genotypes " + refDataType.getName() + ":");
			LOGGER.info("Reference genotypes " + refDataType.getName() + ":");
			for (String path : refBasePaths) {
				System.out.print(" " + path);
				LOGGER.info(" " + path);
			}
			System.out.println();
			System.out.println(" - Reference genotype cache size: " + Ase.DEFAULT_NUMBER_FORMATTER.format(refDataCacheSize));
			LOGGER.info("Reference genotype cache size: " + Ase.DEFAULT_NUMBER_FORMATTER.format(refDataCacheSize));
			if(isSampleToRefSampleFileSet()){
				System.out.println(" - Sample mapping to reference file: " + sampleToRefSampleFile.getAbsolutePath());
				LOGGER.info("Sample mapping to reference file: " + sampleToRefSampleFile.getAbsolutePath());
			}
		}

		if (isGtfSet()) {
			System.out.print(" - GTF file: " + gtf.getAbsolutePath());
			LOGGER.info("GTF file: " + gtf.getAbsolutePath());
		}


		System.out.println();

		System.out.flush(); //flush to make sure config is before errors
		try {
			Thread.sleep(25); //Allows flush to complete
		} catch (InterruptedException ex) {
		}


	}

	public static void printHelp() {
		new HelpFormatter().printHelp(" ", OPTIONS);
	}

	public List<File> getInputFiles() {
		return inputFiles;
	}

	public File getOutputFolder() {
		return outputFolder;
	}

	public int getMinTotalReads() {
		return minTotalReads;
	}

	public int getMinAlleleReads() {
		return minAlleleReads;
	}

	public File getLogFile() {
		return logFile;
	}

	public boolean isDebugMode() {
		return debugMode;
	}

	public int getMinSamples() {
		return minSamples;
	}

	public int getThreads() {
		return threads;
	}

	public String[] getRefBasePaths() {
		return refBasePaths;
	}

	public RandomAccessGenotypeDataReaderFormats getRefDataType() {
		return refDataType;
	}

	public boolean isRefSet() {
		return refBasePaths != null;
	}

	public int getRefDataCacheSize() {
		return refDataCacheSize;
	}

	public File getGtf() {
		return gtf;
	}

	public int getMaxTotalReads() {
		return maxTotalReads;
	}

	public boolean isGtfSet() {
		return gtf != null;
	}

	public double getMinAlleleReadFraction() {
		return minAlleleReadFraction;
	}

	public File getSampleToRefSampleFile() {
		return sampleToRefSampleFile;
	}
	
	public boolean isSampleToRefSampleFileSet(){
		return sampleToRefSampleFile != null;
	}
	
}