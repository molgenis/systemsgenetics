package eqtlmappingpipeline.ase;

import java.io.File;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.ResourceBundle;
import org.apache.commons.cli.ParseException;
import org.apache.log4j.FileAppender;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.SimpleLayout;
import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.multipart.IncompatibleMultiPartGenotypeDataException;
import org.molgenis.genotype.multipart.MultiPartGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.vcf.VcfGenotypeData;

/**
 *
 * @author Patrick Deelen
 */
public class Ase {

	private static final String VERSION = ResourceBundle.getBundle("verion").getString("application.version");
	private static final String HEADER =
			"  /---------------------------------------\\\n"
			+ "  |  Allele Specific Expression Mapper    |\n"
			+ "  |                                       |\n"
			+ "  |             Patrick Deelen            |\n"
			+ "  |        patrickdeelen@gmail.com        |\n"
			+ "  |                                       |\n"
			+ "  | Dasha Zhernakova, Marijke v/d Sijde,  |\n"
			+ "  |   Marc Jan Bonder, Harm-Jan Westra,   |\n"
			+ "  |      Lude Franke, Morris Swertz       |\n"
			+ "  |                                       |\n"
			+ "  |     Genomics Coordication Center      |\n"
			+ "  |        Department of Genetics         |\n"
			+ "  |  University Medical Center Groningen  |\n"
			+ "  \\---------------------------------------/";
	private static final Logger LOGGER = Logger.getLogger(Ase.class);
	private static final DateFormat DATE_TIME_FORMAT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
	private static final Date currentDataTime = new Date();

	public static void main(String[] args) {

		System.out.println(HEADER);
		System.out.println();
		System.out.println("          --- Version: " + VERSION + " ---");
		System.out.println();
		System.out.println("More information: http://molgenis.org/systemsgenetics");
		System.out.println();

		System.out.println("Current date and time: " + DATE_TIME_FORMAT.format(currentDataTime));
		System.out.println();

		System.out.flush(); //flush to make sure header is before errors
		try {
			Thread.sleep(25); //Allows flush to complete
		} catch (InterruptedException ex) {
		}

		if (args.length == 0) {
			AseConfiguration.printHelp();
			System.exit(1);
		}

		final AseConfiguration configuration;
		try {
			configuration = new AseConfiguration(args);
		} catch (ParseException ex) {
			System.err.println("Invalid command line arguments: ");
			System.err.println(ex.getMessage());
			System.err.println();
			AseConfiguration.printHelp();
			System.exit(1);
			return;
		}

		if (!configuration.getOutputFolder().mkdirs()) {
			System.err.println("Failed to create output folder: " + configuration.getOutputFolder().getAbsolutePath());
			System.exit(1);
		}

		startLogging(configuration.getLogFile(), configuration.isDebugMode());
		configuration.printOptions();

		AseResults aseResuls = new AseResults();

		try {

			for (File inputFile : configuration.getInputFiles()) {

				//Loading genotype files:
				GenotypeData genotypeData;
				if (inputFile.isDirectory()) {
					try {
						genotypeData = MultiPartGenotypeData.createFromVcfFolder(inputFile, 100);
					} catch (IncompatibleMultiPartGenotypeDataException ex) {
						System.err.println("Error reading folder with VCF files: " + ex.getMessage());
						LOGGER.fatal("Error reading folder with VCF files: ", ex);
						System.exit(1);
						return;
					}
				} else {
					genotypeData = new VcfGenotypeData(inputFile, 100);
				}

				//TODO test if VCF contains the read depth field

				for (GeneticVariant variant : genotypeData) {

					//Here we are going to do the ASE part

					//Only if variant contains read depth field
					
					//For all samples in file
						
					int a1Count = 0;
					int a2Count = 0;

					if (a1Count + a2Count > configuration.getMinTotalReads()
							&& (a1Count >= configuration.getMinAlleleReads()
							|| a2Count >= configuration.getMinAlleleReads())) {
						
								aseResuls.addResult(variant.getSequenceName(), variant.getStartPos(), variant.getVariantId(), variant.getVariantAlleles().get(0), variant.getVariantAlleles().get(1), a1Count, a2Count);
						
					}


				}


			}

		} catch (IOException ex) {
			System.err.println("Error reading input data: " + ex.getMessage());
			LOGGER.fatal("Error reading input data", ex);
			System.exit(1);
			return;
		} catch (GenotypeDataException ex) {
			System.err.println("Error reading input data: " + ex.getMessage());
			LOGGER.fatal("Error reading input data", ex);
			System.exit(1);
			return;
		}
		
		for(AseVariant aseVariant : aseResuls ){
			
			
			
		}
		

	}

	private static void startLogging(File logFile, boolean debugMode) {
		try {
			FileAppender logAppender = new FileAppender(new SimpleLayout(), logFile.getCanonicalPath(), false);
			Logger.getRootLogger().removeAllAppenders();
			Logger.getRootLogger().addAppender(logAppender);
			if (debugMode) {
				LOGGER.setLevel(Level.DEBUG);
			} else {
				LOGGER.setLevel(Level.INFO);
			}
		} catch (IOException e) {
			System.err.println("Failed to create logger: " + e.getMessage());
			System.exit(1);
		}



		LOGGER.info(
				"\n" + HEADER);
		LOGGER.info("Version: " + VERSION);
		LOGGER.info("Current date and time: " + DATE_TIME_FORMAT.format(currentDataTime));
		LOGGER.info("Log level: " + LOGGER.getLevel());

		System.out.println("Started logging");
		System.out.println();
	}
}
