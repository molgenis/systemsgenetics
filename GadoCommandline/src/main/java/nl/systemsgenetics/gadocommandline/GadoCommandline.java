/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.gadocommandline;

import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.ResourceBundle;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.Appender;
import org.apache.logging.log4j.core.LoggerContext;
import org.apache.logging.log4j.core.appender.ConsoleAppender;
import org.apache.logging.log4j.core.appender.FileAppender;
import org.apache.logging.log4j.core.config.Configuration;
import org.apache.logging.log4j.core.config.Configurator;
import org.apache.logging.log4j.core.config.LoggerConfig;
import org.apache.logging.log4j.core.layout.PatternLayout;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.multipart.IncompatibleMultiPartGenotypeDataException;

/**
 *
 * @author patri
 */
public class GadoCommandline {

	public static final String VERSION = ResourceBundle.getBundle("verion").getString("application.version");
	private static final DateFormat DATE_TIME_FORMAT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
	private static final Logger LOGGER = LogManager.getLogger(GadoCommandline.class);
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
		} catch (InterruptedException ex) {
		}

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

			if (options.isDebugMode()) {
				Configurator.setRootLevel(Level.DEBUG);
			} else {
				Configurator.setRootLevel(Level.INFO);
			}
			LoggerContext context = LoggerContext.getContext(false);
			Configuration config = context.getConfiguration();

			PatternLayout loggingLayoutFull = PatternLayout.newBuilder()
					.withPattern("[%level] %d{ABSOLUTE} - %c{1} - %msg%n")
					.build();

			PatternLayout loggingLayoutReduced = PatternLayout.newBuilder()
					.withPattern("%msg%n")
					.build();

			// Stdout appender
			ConsoleAppender stdOut = ConsoleAppender.newBuilder()
					.setName("stdout")
					.setLayout(loggingLayoutReduced)
					.build();
			stdOut.start();

			// Log file appender
			FileAppender file = FileAppender.newBuilder()
					.setName("file")
					.setLayout(loggingLayoutFull)
					.withFileName(options.getLogFile().getCanonicalPath())
					.build();
			file.start();

			// Make sure any existing loggers are removed
			for (Appender appender : context.getRootLogger().getAppenders().values()) {
				context.getRootLogger().removeAppender(appender);
			}

			// Add the appenders to the root logger
			Logger rootLogger = context.getRootLogger();
			LoggerConfig rootLoggerConfig = config.getRootLogger();

			rootLoggerConfig.addAppender(stdOut, Level.INFO, null);
			rootLoggerConfig.addAppender(file, options.isDebugMode() ? Level.DEBUG : Level.INFO, null);
			config.addLogger(rootLogger.getName(), rootLoggerConfig);

			context.updateLoggers(config);
			
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
				case EXPAND_PREDICTIONS:
					SpikeInKnownAnnotations.spikeIn(options);
					break;
				case PREPARE_HPO_FOR_PREDICTIONS:
					PrepareHpoPhenoToGenesForPredictions.prioritize(options);
					break;
				case CONVERT_TXT:
					Converters.convertTxtToBin(options);;
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
