/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.ase;

import java.io.File;
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

/**
 *
 * @author Patrick Deelen
 */
public class AseConfiguration {

	private static final Logger LOGGER;
	private static final Options OPTIONS;
	private final List<File> inputFiles;
	private final File outputFolder;
	private final int minTotalReads;
	private final int minAlleleReads;
	private final File logFile;
	private final boolean debugMode;

	static {

		LOGGER = Logger.getLogger(AseConfiguration.class);

		OPTIONS = new Options();

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Path to a vcf.gz file or a folder with per chr 1 vcf.gz file. Tabix file must be present");
		OptionBuilder.withLongOpt("input");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create('i'));

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
		OptionBuilder.withDescription("Min number of reads per alle");
		OptionBuilder.withLongOpt("minAlleleReads");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create('a'));

		OptionBuilder.withArgName("boolean");
		OptionBuilder.withDescription("Activate debug mode. This will result in a more verbose log file");
		OptionBuilder.withLongOpt("debug");
		OPTIONS.addOption(OptionBuilder.create('d'));

	}

	public AseConfiguration(String... args) throws ParseException {

		final CommandLine commandLine = new PosixParser().parse(OPTIONS, args, true);

		String[] inputPaths = commandLine.getOptionValues('i');
		ArrayList<File> inputFilesTmp = new ArrayList<File>(inputPaths.length);

		for (String inputPath : inputPaths) {
			inputFilesTmp.add(new File(inputPath));
		}

		inputFiles = Collections.unmodifiableList(inputFilesTmp);

		outputFolder = new File(commandLine.getOptionValue('o'));

		logFile = new File(outputFolder, "ase.log");

		try {
			minTotalReads = Integer.parseInt(commandLine.getOptionValue('r'));
		} catch (NumberFormatException e) {
			throw new ParseException("Error parsing --minReads \"" + commandLine.getOptionValue('r') + "\" is not an int");
		}

		try {
			minAlleleReads = Integer.parseInt(commandLine.getOptionValue('a'));
		} catch (NumberFormatException e) {
			throw new ParseException("Error parsing --minAlleleReads \"" + commandLine.getOptionValue('a') + "\" is not an int");
		}
		
		debugMode = commandLine.hasOption('d');


	}
	
	public void printOptions() {

		System.out.println("Interpreted arguments: ");
		System.out.println(" - Input base path: ");
		LOGGER.info("Input base path: ");
		
		for(File inputFile : inputFiles){
			System.out.println("  * " + inputFile.getAbsolutePath());
			LOGGER.info(" * " + inputFile.getAbsolutePath());
		}
		
		System.out.println(" - Output folder: " + outputFolder.getAbsolutePath());
		LOGGER.info("Output folder: " + outputFolder.getAbsolutePath());
		
		System.out.println(" - Minimum number of reads per genotype: " + minTotalReads);
		LOGGER.info("Minimum number of reads per genotype: " + minTotalReads);		
		
		System.out.println(" - Minimum number of reads per allele: " + minAlleleReads);
		LOGGER.info("Minimum number of reads per allele: " + minAlleleReads);
		
		LOGGER.debug("Debug mode activated");


		System.out.println();


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
	
	
}
