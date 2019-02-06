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
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;

/**
 *
 * @author patri
 */
public class Depict2Options {

	private static final Options OPTIONS;
	private static int numberOfThreadsToUse = Runtime.getRuntime().availableProcessors();

	private final String[] genotypeBasePath;
	private final RandomAccessGenotypeDataReaderFormats genotypeType;
	private final File outputFile;
	private final String gwasZscoreMatrixPath;
	private final int numberOfPermutations;
	private final int windowExtend;
	private final double maxRBetweenVariants;

	static {

		OPTIONS = new Options();

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("GWAS Z-sccore binary matrix. Rows variants, Cols phenotypes. Without .dat");
		OptionBuilder.withLongOpt("gwas");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("g"));

		OptionBuilder.withArgName("basePath");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("The reference genotypes");
		OptionBuilder.withLongOpt("referenceGenotypes");
		OPTIONS.addOption(OptionBuilder.create("r"));

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
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("p"));

		OptionBuilder.withArgName("int");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Number of bases to add left and right of gene window");
		OptionBuilder.withLongOpt("window");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("w"));

		OptionBuilder.withArgName("double");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Max correlation between variants to use (recommend = 0.95)");
		OptionBuilder.withLongOpt("variantCorrelation");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("v"));

	}

	public Depict2Options(String... args) throws ParseException {

		CommandLineParser parser = new PosixParser();
		final CommandLine commandLine = parser.parse(OPTIONS, args, false);

//		eigenVectorFile = new File(commandLine.getOptionValue("e"));
//		pathwayMatrixFile = new File(commandLine.getOptionValue("p"));
//
//		predictionsFile = new File(commandLine.getOptionValue("o", eigenVectorFile.getPath() + ".GenesetZScores.txt"));
//
//		backgroudGenesFile = commandLine.hasOption("b") ? new File(commandLine.getOptionValue("b")) : null;
		if (commandLine.hasOption('t')) {
			try {
				numberOfThreadsToUse = Integer.parseInt(commandLine.getOptionValue('t'));
			} catch (NumberFormatException e) {
				throw new ParseException("Error parsing --threads \"" + commandLine.getOptionValue('t') + "\" is not an int");
			}
		}

		try {
			numberOfPermutations = Integer.parseInt(commandLine.getOptionValue('p'));
		} catch (NumberFormatException e) {
			throw new ParseException("Error parsing --permutations \"" + commandLine.getOptionValue('p') + "\" is not an int");
		}

		try {
			windowExtend = Integer.parseInt(commandLine.getOptionValue('w'));
		} catch (NumberFormatException e) {
			throw new ParseException("Error parsing --window \"" + commandLine.getOptionValue('w') + "\" is not an int");
		}

		try {
			maxRBetweenVariants = Double.parseDouble(commandLine.getOptionValue('v'));
		} catch (NumberFormatException e) {
			throw new ParseException("Error parsing --variantCorrelation \"" + commandLine.getOptionValue('v') + "\" is not an double");
		}

		outputFile = new File(commandLine.getOptionValue('o'));
		gwasZscoreMatrixPath = commandLine.getOptionValue('g');

		if (commandLine.hasOption('r')) {

			genotypeBasePath = commandLine.getOptionValues('r');

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

		} else {
			genotypeBasePath = null;
			genotypeType = null;
		}

	}

	public static void printHelp() {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp(" ", OPTIONS);
	}

	public void printOptions() {

		StringBuilder genotypeBasePaths = new StringBuilder();
		for (String path : genotypeBasePath) {
			genotypeBasePaths.append(path);
			genotypeBasePaths.append(' ');
		}

		if (genotypeBasePath != null) {
			System.out.println(" - Reference genotype data: " + genotypeBasePaths);
			System.out.println(" - Reference genotype data type: " + genotypeType.getName());
		}

		System.out.println(" - Ouput path: " + outputFile.getPath());
		
		System.out.println(" - Gwas Z-score matrix: " + gwasZscoreMatrixPath);
		
		System.out.println(" - Gene window extend in bases: " + windowExtend);
		System.out.println(" - Number of permutations: " + numberOfPermutations);
		System.out.println(" - Max correlation between variants: " + maxRBetweenVariants);
		System.out.println(" - Number of threads to use: " + numberOfThreadsToUse);
		
		
		System.out.println(" - ");
		
		
	}

	public String[] getGenotypeBasePath() {
		return genotypeBasePath;
	}

	public RandomAccessGenotypeDataReaderFormats getGenotypeType() {
		return genotypeType;
	}

	public File getOutputFile() {
		return outputFile;
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

}
