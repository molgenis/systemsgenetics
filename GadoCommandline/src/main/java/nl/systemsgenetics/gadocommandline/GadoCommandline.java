/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.gadocommandline;

import java.io.File;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.ResourceBundle;
import org.apache.commons.cli.ParseException;
import org.apache.log4j.ConsoleAppender;
import org.apache.log4j.FileAppender;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.SimpleLayout;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.multipart.IncompatibleMultiPartGenotypeDataException;
import org.molgenis.genotype.tabix.TabixFileNotFoundException;

/**
 *
 * @author patri
 */
public class GadoCommandline {

	public static final String VERSION = ResourceBundle.getBundle("verion").getString("application.version");
	private static final DateFormat DATE_TIME_FORMAT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
	private static final Logger LOGGER = Logger.getLogger(GadoCommandline.class);
	private static final String HEADER
			= "  /---------------------------------------\\\n"
			+ "  |                 GADO                  |\n"
			+ "  |                                       |\n"
			+ "  |  University Medical Center Groningen  |\n"
			+ "  \\---------------------------------------/";

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) {
		
		System.out.println(HEADER);
		System.out.println();
		System.out.println("          --- Version: " + VERSION + " ---");
		System.out.println();
		System.out.println("More information: http://molgenis.org/systemsgenetics");
		System.out.println();

		Date currentDataTime = new Date();
		String startDateTime = DATE_TIME_FORMAT.format(currentDataTime);
		System.out.println("Current date and time: " + startDateTime);
		System.out.println();

		System.out.flush(); //flush to make sure header is before errors
		try {
			Thread.sleep(25); //Allows flush to complete
		} catch (InterruptedException ex) {}

		GadoOptions options;

		if (args.length == 0) {
			GadoOptions.printHelp();
			return;
		}

		try {
			options = new GadoOptions(args);
		} catch (ParseException ex) {
			System.err.println("Error parsing commandline: " + ex.getMessage());
			GadoOptions.printHelp();
			return;
		}

		if (options.getLogFile().getParentFile() != null && !options.getLogFile().getParentFile().isDirectory()) {
			if (!options.getLogFile().getParentFile().mkdirs()) {
				System.err.println("Failed to create output folder: " + options.getLogFile().getParent());
				System.exit(1);
			}
		}
		
		try {
			FileAppender logFileAppender = new FileAppender(new SimpleLayout(), options.getLogFile().getCanonicalPath(), false);
			ConsoleAppender logConsoleInfoAppender = new ConsoleAppender(new InfoOnlyLogLayout());
			Logger.getRootLogger().removeAllAppenders();
			Logger.getRootLogger().addAppender(logFileAppender);

			LOGGER.info("DEPICT" + VERSION);
			LOGGER.info("Current date and time: " + startDateTime);

			Logger.getRootLogger().addAppender(logConsoleInfoAppender);

			if (options.isDebugMode()) {
				Logger.getRootLogger().setLevel(Level.DEBUG);
			} else {
				Logger.getRootLogger().setLevel(Level.INFO);
			}

		} catch (IOException e) {
			System.err.println("Failed to create logger: " + e.getMessage());
			System.exit(1);
		}

		options.printOptions();

		try {
			switch (options.getMode()) {
				case PROCESS:
					ProcessCaseHpo.process(options);
					break;
				case PRIORITIZE:
					HpoGenePrioritisation.prioritize(options);
					break;
			}
		} catch (IOException e) {
			System.err.println("Problem running mode: " + options.getMode());
			System.err.println("Error accessing input data: " + e.getMessage());
			System.err.println("See log file for stack trace");
			LOGGER.fatal("Error accessing input data: " + e.getMessage(), e);
			System.exit(1);
		} catch (IncompatibleMultiPartGenotypeDataException e) {
			System.err.println("Problem running mode: " + options.getMode());
			System.err.println("Error combining the impute genotype data files: " + e.getMessage());
			System.err.println("See log file for stack trace");
			LOGGER.fatal("Error combining the impute genotype data files: " + e.getMessage(), e);
			System.exit(1);
		} catch (GenotypeDataException e) {
			System.err.println("Problem running mode: " + options.getMode());
			System.err.println("Error reading genotype data: " + e.getMessage());
			System.err.println("See log file for stack trace");
			LOGGER.fatal("Error reading input data: " + e.getMessage(), e);
			System.exit(1);
		} catch (Exception e) {
			System.err.println("Problem running mode: " + options.getMode());
			System.err.println("Error message: " + e.getMessage());
			System.err.println("See log file for stack trace");
			LOGGER.fatal("Error: " + e.getMessage(), e);
			if (LOGGER.isDebugEnabled()) {
				e.printStackTrace();
			}
			System.exit(1);
		}
		
		LOGGER.info("Analysis completed");
		
	}

}
