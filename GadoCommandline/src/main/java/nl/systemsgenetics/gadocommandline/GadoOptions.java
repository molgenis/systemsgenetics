/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.gadocommandline;

import java.io.File;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;

/**
 *
 * @author patri
 */
public class GadoOptions {

	private static final Options OPTIONS;
	private static final Logger LOGGER = LogManager.getLogger(GadoOptions.class);

	private final GadoMode mode;

	private final File logFile;
	private final boolean debugMode;

	private final File outputBasePath;
	private final File genesFile;
	private final File hpoOboFile;
	private final File predictionMatrixFile;
	private final File predictionInfoFile;
	private final File inputCaseHpoFile;
	private final File processedCaseHpoFile;
	private final File annotationMatrixFile;
	private final File phenotypesToGenesFile;

	public boolean isDebugMode() {
		return debugMode;
	}

	static {

		OPTIONS = new Options();

		OptionBuilder.withArgName("string");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("One of the following modes:\n"
				+ "* PROCESS - Process the HPO terms of cases. Suggests  parent terms if needed.\n"
				+ "* PRIORITIZE - Uses output of PROCESS to prioritize genes\n"
				+ "* EXPAND_PREDICTIONS - Spike-in known HPO-Gene annotation into the prediction matrix. \n"
				+ "* PREPARE_HPO_FOR_PREDICTIONS - Converts the HPO file \"phenotype_to_genes.txt\" to a matrix that can be used in gene network predictions");
		OptionBuilder.withLongOpt("mode");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("m"));

		OptionBuilder.withArgName("file");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("HPO terms per case. Single line per case. First col is case ID, followed by tab separated  HPO terms.");
		OptionBuilder.withLongOpt("caseHpo");
		OPTIONS.addOption(OptionBuilder.create("ch"));

		OptionBuilder.withArgName("file");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("The output of mode PROCESS. Type x in the last column to exclude a term.");
		OptionBuilder.withLongOpt("caseHpoProcessed");
		OPTIONS.addOption(OptionBuilder.create("chp"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("The output path. For mode PRIORITIZE supply a folder for the output files per case. For mode EXPAND_PREDICTIONS output matrix");
		OptionBuilder.withLongOpt("output");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("o"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("File with gene info. col1: geneName (ensg) col2: HGNC symbol");
		OptionBuilder.withLongOpt("genes");
		OPTIONS.addOption(OptionBuilder.create("g"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("HPO ontology file");
		OptionBuilder.withLongOpt("hpoOntology");
		OPTIONS.addOption(OptionBuilder.create("ho"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("HPO prediction matrix in binary format. In the case of mode CONVERT_TXT this is used to point to the .txt matrix file.");
		OptionBuilder.withLongOpt("hpoPredictions");
		OPTIONS.addOption(OptionBuilder.create("hp"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("HPO annotion matrix in txt or txt.gz format. Only for mode EXPAND_PREDICTIONS Rows: genes, cols: hpo. Must be exact same row and col order as prediction matrix ");
		OptionBuilder.withLongOpt("hpoAnnotations");
		OPTIONS.addOption(OptionBuilder.create("hpa"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("HPO predictions info");
		OptionBuilder.withLongOpt("hpoPredictionsInfo");
		OPTIONS.addOption(OptionBuilder.create("hpi"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("HPO predictions info");
		OptionBuilder.withLongOpt("phenotypesToGenes");
		OPTIONS.addOption(OptionBuilder.create("ptg"));

//		OptionBuilder.withArgName("boolean");
//		OptionBuilder.withDescription("Activate debug mode. This will result in a more verbose log file and will save many intermediate results to files. Not recommended for large analysis.");
//		OptionBuilder.withLongOpt("debug");
//		OPTIONS.addOption(OptionBuilder.create("d"));
	}

	public GadoOptions(String... args) throws ParseException {

		final CommandLineParser parser = new PosixParser();
		final CommandLine commandLine = parser.parse(OPTIONS, args, false);

		outputBasePath = new File(commandLine.getOptionValue('o'));
		logFile = new File(outputBasePath + "gado.log");
		debugMode = false; //commandLine.hasOption('d');

		try {
			mode = GadoMode.valueOf(commandLine.getOptionValue("m").toUpperCase());
		} catch (IllegalArgumentException e) {
			throw new ParseException("Error parsing --mode \"" + commandLine.getOptionValue("m") + "\" is not a valid mode");
		}

		if (commandLine.hasOption("hpi")) {
			predictionInfoFile = new File(commandLine.getOptionValue("hpi"));
		} else {
			if (mode == GadoMode.PROCESS) {
				throw new ParseException("Option --hpoPredictionsInfo is requered for mode PROCESS");
			}
			predictionInfoFile = null;
		}

		if (commandLine.hasOption("ho")) {
			hpoOboFile = new File(commandLine.getOptionValue("ho"));
		} else {
			if (mode == GadoMode.PROCESS) {
				throw new ParseException("Option --hpoOntology is requered for mode PROCESS");
			}
			hpoOboFile = null;
		}

		if (commandLine.hasOption("ch")) {
			inputCaseHpoFile = new File(commandLine.getOptionValue("ch"));
		} else {
			if (mode == GadoMode.PROCESS) {
				throw new ParseException("Option --caseHpo is requered for mode PROCESS");
			}
			inputCaseHpoFile = null;
		}

		if (commandLine.hasOption("g")) {
			genesFile = new File(commandLine.getOptionValue("g"));
		} else {
			if (mode == GadoMode.PRIORITIZE || mode == GadoMode.PREPARE_HPO_FOR_PREDICTIONS) {
				throw new ParseException("Option --genes is requered for mode PRIORITIZE and PREPARE_HPO_FOR_PREDICTIONS");
			}
			genesFile = null;
		}

		if (commandLine.hasOption("hp")) {
			String hp = commandLine.getOptionValue("hp");
			if (hp.endsWith(".dat")) {
				hp = hp.substring(0, hp.length() - 4);
			}
			predictionMatrixFile = new File(hp);
		} else {
			if (mode == GadoMode.PRIORITIZE) {
				throw new ParseException("Option --hpoPredictions is requered for mode PRIORITIZE");
			}
			if (mode == GadoMode.CONVERT_TXT) {
				throw new ParseException("Option --hpoPredictions is requered for mode CONVERT_TXT");
			}
			predictionMatrixFile = null;
		}

		if (commandLine.hasOption("hpa")) {
			String hp = commandLine.getOptionValue("hpa");
			if (hp.endsWith(".dat")) {
				hp = hp.substring(0, hp.length() - 4);
			}
			annotationMatrixFile = new File(hp);
		} else {
			if (mode == GadoMode.EXPAND_PREDICTIONS) {
				throw new ParseException("Option --hpoAnnotations is requered for mode EXPAND_PREDICTIONS");
			}
			annotationMatrixFile = null;
		}

		if (commandLine.hasOption("chp")) {
			processedCaseHpoFile = new File(commandLine.getOptionValue("chp"));
		} else {
			if (mode == GadoMode.PRIORITIZE) {
				throw new ParseException("Option --caseHpoProcessed is requered for mode PRIORITIZE");
			}
			processedCaseHpoFile = null;
		}
		
		if (commandLine.hasOption("ptg")) {
			phenotypesToGenesFile = new File(commandLine.getOptionValue("ptg"));
		} else {
			if (mode == GadoMode.PREPARE_HPO_FOR_PREDICTIONS) {
				throw new ParseException("Option --phenotypesToGenes is requered for mode PREPARE_HPO_FOR_PREDICTIONS");
			}
			phenotypesToGenesFile = null;
		}

	}

	public static void printHelp() {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp(" ", OPTIONS);
	}

	public void printOptions() {

		LOGGER.info("Supplied options:");

		LOGGER.info(" * Mode: " + mode.name());
		LOGGER.info(" * Output: " + outputBasePath);
		switch (mode) {
			case PROCESS:
				LOGGER.info(" * Case HPO terms: " + this.inputCaseHpoFile.getAbsolutePath());
				LOGGER.info(" * HPO ontology file: " + this.hpoOboFile.getAbsolutePath());
				LOGGER.info(" * Prediction info file: " + this.predictionInfoFile.getAbsolutePath());
				break;
			case PRIORITIZE:
				LOGGER.info(" * Processed case HPO terms:" + this.processedCaseHpoFile.getAbsolutePath());
				LOGGER.info(" * HPO prediction matrix: " + this.predictionMatrixFile.getAbsolutePath());
				LOGGER.info(" * Genes info: " + this.genesFile.getAbsolutePath());
				break;
			case EXPAND_PREDICTIONS:
				LOGGER.info(" * HPO prediction matrix: " + this.predictionMatrixFile.getAbsolutePath());
				LOGGER.info(" * HPO annotation matrix: " + this.annotationMatrixFile.getAbsolutePath());
				break;
			case PREPARE_HPO_FOR_PREDICTIONS:
				LOGGER.info(" * Phenotype to genes files: " + this.phenotypesToGenesFile.getAbsolutePath());
				LOGGER.info(" * Genes info: " + this.genesFile.getAbsolutePath());
				break;
			case CONVERT_TXT:
				LOGGER.info(" * HPO prediction matrix in txt format: " + this.predictionMatrixFile.getAbsolutePath());
				break;
		}

		//LOGGER.info(" * Debug mode: " + (debugMode ? "on" : "off"));
	}

	public File getLogFile() {
		return logFile;
	}

	public GadoMode getMode() {
		return mode;
	}

	public String getOutputBasePath() {
		return outputBasePath.getPath();
	}

	public File getGenesFile() {
		return genesFile;
	}

	public File getHpoOboFile() {
		return hpoOboFile;
	}

	public File getPredictionMatrixFile() {
		return predictionMatrixFile;
	}

	public File getPredictionInfoFile() {
		return predictionInfoFile;
	}

	public File getInputCaseHpoFile() {
		return inputCaseHpoFile;
	}

	public File getProcessedCaseHpoFile() {
		return processedCaseHpoFile;
	}

	public File getAnnotationMatrixFile() {
		return annotationMatrixFile;
	}

	public File getPhenotypeToGenes() {
		return phenotypesToGenesFile;
	}

}
