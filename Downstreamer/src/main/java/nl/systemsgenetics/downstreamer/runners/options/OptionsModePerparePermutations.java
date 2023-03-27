/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.runners.options;

import java.io.File;
import static nl.systemsgenetics.downstreamer.runners.options.OptionsBase.OPTIONS;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 *
 * @author patri
 */
public class OptionsModePerparePermutations extends OptionsBase {

	private static final Logger LOGGER = LogManager.getLogger(OptionsBase.class);

	private final boolean forceNormalGenePvalues;
	private final boolean special;
	private final File geneInfoFile;
	private final File permutationFolder;
	private final File chunkFile;

	static {

		OptionBuilder.withArgName("boolean");
		OptionBuilder.withDescription("Force normal gene p-values before pathway enrichtment");
		OptionBuilder.withLongOpt("forceNormalGenePfngpvalues");
		OPTIONS.addOption(OptionBuilder.create("fngp"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("File with gene info. col1: geneName (ensg) col2: chr col3: startPos col4: stopPos col5: geneType col6: chrArm");
		OptionBuilder.withLongOpt("genes");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("ge"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Folder with permuataion results");
		OptionBuilder.withLongOpt("permutations");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("p"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("File with chunks");
		OptionBuilder.withLongOpt("chunks");
		OPTIONS.addOption(OptionBuilder.create("c"));

		OptionBuilder.withArgName("boolean");
		OptionBuilder.withDescription("Special");
		OptionBuilder.withLongOpt("special");
		OPTIONS.addOption(OptionBuilder.create("s"));

	}

	public OptionsModePerparePermutations(String[] args) throws ParseException {
		super(args);

		// Parse arguments
		final CommandLineParser parser = new PosixParser();
		final CommandLine commandLine = parser.parse(OPTIONS, args, false);

		forceNormalGenePvalues = commandLine.hasOption("fngp");
		special = commandLine.hasOption("s");
		geneInfoFile = new File(commandLine.getOptionValue("ge"));
		permutationFolder = new File(commandLine.getOptionValue("p"));
		if(!special){
			chunkFile = new File(commandLine.getOptionValue("c"));
		} else {
			chunkFile = null;
		}
		

		printOptions();

	}

	@Override
	public void printOptions() {
		super.printOptions();

		LOGGER.info(" * geneInfoFile: " + geneInfoFile.getPath());
		if (!special) {
			LOGGER.info(" * chunk file: " + chunkFile.getPath());
		} else {
			LOGGER.info("Special mode");
		}
		LOGGER.info(" * permutation folder: " + permutationFolder.getPath());
		LOGGER.info(" * forceNormalGenePvalues: " + forceNormalGenePvalues);
		

	}

	public boolean isForceNormalGenePvalues() {
		return forceNormalGenePvalues;
	}

	public File getGeneInfoFile() {
		return geneInfoFile;
	}

	public File getPermutationFolder() {
		return permutationFolder;
	}

	public File getChunkFile() {
		return chunkFile;
	}

	public boolean isSpecial() {
		return special;
	}
}
