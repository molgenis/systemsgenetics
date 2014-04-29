package eqtlmappingpipeline.ase;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;
import java.util.Iterator;
import java.util.List;
import org.apache.commons.cli.ParseException;
import org.apache.log4j.FileAppender;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.SimpleLayout;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.multipart.IncompatibleMultiPartGenotypeDataException;
import org.molgenis.genotype.multipart.MultiPartGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.GeneticVariantMeta;
import org.molgenis.genotype.variant.GenotypeRecord;
import org.molgenis.genotype.vcf.VcfGenotypeData;

/**
 *
 * @author Patrick Deelen
 */
public class Ase {

//	private static final String VERSION = ResourceBundle.getBundle("verion").getString("application.version");
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
//		System.out.println("          --- Version: " + VERSION + " ---");
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

		if (!configuration.getOutputFolder().isDirectory() && !configuration.getOutputFolder().mkdirs()) {
			System.err.println("Failed to create output folder: " + configuration.getOutputFolder().getAbsolutePath());
			System.exit(1);
		}

		startLogging(configuration.getLogFile(), configuration.isDebugMode());
		configuration.printOptions();

		AseResults aseResults = new AseResults();
		int sampleCounter = 0;
		int fileCounter = 0;

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

					if (variant.getVariantMeta().getRecordType("AD") == GeneticVariantMeta.Type.INTEGER_LIST) {
						//include variant

						Iterator<Sample> sampleIterator = genotypeData.getSamples().iterator();

						for (GenotypeRecord record : variant.getSampleGenotypeRecords()) {

							Sample sample = sampleIterator.next();

							if (!record.containsGenotypeRecord("AD")) {
								continue;
							}

							try {
								Alleles alleles = record.getSampleAlleles();
								if(alleles.getAlleleCount() != 2 || alleles.get(0) == alleles.get(1)){
									continue;
								}
								
								List<Integer> counts = (List<Integer>) record.getGenotypeRecordData("AD");

								int a1Count = counts.get(0);
								int a2Count = counts.get(1);

								if (a1Count + a2Count > configuration.getMinTotalReads()
										&& (a1Count >= configuration.getMinAlleleReads()
										|| a2Count >= configuration.getMinAlleleReads())) {

									aseResults.addResult(variant.getSequenceName(), variant.getStartPos(), variant.getVariantId(), variant.getVariantAlleles().get(0), variant.getVariantAlleles().get(1), a1Count, a2Count);

								}

							} catch (GenotypeDataException ex) {
								System.err.println("Error parsing " + variant.getSequenceName() + ":" + variant.getStartPos() + " for sample " + sample.getId() + " " + ex.getMessage());
								LOGGER.fatal("Error parsing " + variant.getSequenceName() + ":" + variant.getStartPos() + " for sample " + sample.getId(), ex);
								System.exit(1);
								return;
							}


						}

					}

				}

				fileCounter += 1;
				
				if(fileCounter % 100 == 0){
					System.out.println("Loaded "  + fileCounter + " out of " + configuration.getInputFiles().size() + " files");
				}
				
				sampleCounter += genotypeData.getSampleNames().length;

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
		
		LOGGER.info("Loading files complete. " + sampleCounter + " samples encountered");
		System.out.println("Loading files complete. " + sampleCounter + " samples encountered");
		
		Iterator<AseVariant> aseIterator = aseResults.iterator();
		while(aseIterator.hasNext()){
			if(aseIterator.next().getSampleCount() < configuration.getMinSamples()){
				aseIterator.remove();
			}
		}
		
		AseVariant[] aseVariants = new AseVariant[aseResults.getCount()];
		{
			int i = 0;
			for (AseVariant aseVariant : aseResults) {

				//This can be made multithreaded if needed
				aseVariant.calculateMetaZscoreAndPvalue();

				aseVariants[i] = aseVariant;
				++i;

			}
		}


		Arrays.sort(aseVariants);

		File outputFile = new File(configuration.getOutputFolder(), "result.txt");
		try {

			printAseResults(outputFile, aseVariants);

		} catch (UnsupportedEncodingException ex) {
			throw new RuntimeException(ex);
		} catch (FileNotFoundException ex) {
			System.err.println("Unable to create output file at " + outputFile.getAbsolutePath());
			LOGGER.fatal("Unable to create output file at " + outputFile.getAbsolutePath(), ex);
			System.exit(1);
			return;
		} catch (IOException ex) {
			System.err.println("Unable to create output file at " + outputFile.getAbsolutePath());
			LOGGER.fatal("Unable to create output file at " + outputFile.getAbsolutePath(), ex);
			System.exit(1);
			return;
		}

		System.out.println("Program completed");
		LOGGER.info("Program completed");


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
//		LOGGER.info("Version: " + VERSION);
		LOGGER.info("Current date and time: " + DATE_TIME_FORMAT.format(currentDataTime));
		LOGGER.info("Log level: " + LOGGER.getLevel());

		System.out.println("Started logging");
		System.out.println();
	}

	/**
	 * 
	 * @param outputFile
	 * @param aseVariants
	 * @throws UnsupportedEncodingException
	 * @throws FileNotFoundException
	 * @throws IOException 
	 */
	private static void printAseResults(File outputFile, AseVariant[] aseVariants) throws UnsupportedEncodingException, FileNotFoundException, IOException {
		
		BufferedWriter outputWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outputFile), AseConfiguration.ENCODING));

		outputWriter.append("Meta_P\tMeta_Z\tChr\tPos\tSnpId\tSample_Count\tRef_Allele\tAlt_Allele\tRef_Counts\tAlt_Counts\n");


		for (AseVariant aseVariant : aseVariants) {

			outputWriter.append(String.valueOf(aseVariant.getMetaPvalue()));
			outputWriter.append('\t');
			outputWriter.append(String.valueOf(aseVariant.getMetaZscore()));
			outputWriter.append('\t');
			outputWriter.append(aseVariant.getChr());
			outputWriter.append('\t');
			outputWriter.append(String.valueOf(aseVariant.getPos()));
			outputWriter.append('\t');
			outputWriter.append(String.valueOf(aseVariant.getSampleCount()));
			outputWriter.append('\t');
			outputWriter.append(aseVariant.getId().getPrimairyId() == null ? "." : aseVariant.getId().getPrimairyId());
			outputWriter.append('\t');
			outputWriter.append(aseVariant.getA1().getAlleleAsString());
			outputWriter.append('\t');
			outputWriter.append(aseVariant.getA2().getAlleleAsString());
			outputWriter.append('\t');

			for (int i = 0; i < aseVariant.getA1Counts().size(); ++i) {
				if (i > 0) {
					outputWriter.append(',');
				}
				outputWriter.append(String.valueOf(aseVariant.getA1Counts().getQuick(i)));
			}
			outputWriter.append('\t');
			for (int i = 0; i < aseVariant.getA2Counts().size(); ++i) {
				if (i > 0) {
					outputWriter.append(',');
				}
				outputWriter.append(String.valueOf(aseVariant.getA2Counts().getQuick(i)));
			}
			outputWriter.append('\n');

		}


		outputWriter.close();
	}
}
