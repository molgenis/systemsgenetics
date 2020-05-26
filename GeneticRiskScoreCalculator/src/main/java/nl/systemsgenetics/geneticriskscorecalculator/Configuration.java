package nl.systemsgenetics.geneticriskscorecalculator;

import java.io.File;
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
public class Configuration {

	private static final Logger LOGGER;
	private static final Options OPTIONS;
	private final String[] inputPaths;
	private final RandomAccessGenotypeDataReaderFormats inputDataType;
	private final File outputFile;
	private final File risksnpsFile;
        private final String onlyCount;
        private final String harmonizedData;
        private final double inclusionThreshold;
        
	static {

		LOGGER = Logger.getLogger(Configuration.class);

		OPTIONS = new Options();

                // isRequired gebruiken ja of nee
                // commandline.hash option om te kijken of ie bestaat
                
                
		OptionBuilder.withArgName("basePath");
		OptionBuilder.hasArgs();
		OptionBuilder.isRequired();
		OptionBuilder.withDescription("The path to the genotype data used for the risk score calculation");
		OptionBuilder.withLongOpt("input");
		OPTIONS.addOption(OptionBuilder.create("i"));

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
		OptionBuilder.withLongOpt("inputType");
		OPTIONS.addOption(OptionBuilder.create("I"));

                OptionBuilder.withArgName("risksnps");
		OptionBuilder.hasArgs();
		OptionBuilder.isRequired();
		OptionBuilder.withDescription("The input file with risk snps");
		OptionBuilder.withLongOpt("risksnps");
		OPTIONS.addOption(OptionBuilder.create("r"));
                
		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.isRequired();
		OptionBuilder.withDescription("The output file with risk scores");
		OptionBuilder.withLongOpt("output");
		OPTIONS.addOption(OptionBuilder.create("o"));

		OptionBuilder.withArgName("count");
		OptionBuilder.hasArgs();                
		OptionBuilder.withDescription("Only count risk alleles");
		OptionBuilder.withLongOpt("count");
		OPTIONS.addOption(OptionBuilder.create("c"));
                
		OptionBuilder.withArgName("threshold");
		OptionBuilder.hasArgs();                
		OptionBuilder.isRequired();
                OptionBuilder.withDescription("P value threshold for SNP inclusion");
		OptionBuilder.withLongOpt("thres");
		OPTIONS.addOption(OptionBuilder.create("t"));

		OptionBuilder.withArgName("harmonized");
		OptionBuilder.hasArgs();                
                OptionBuilder.withDescription("Genotype data and risk alleles are on same strand");
		OptionBuilder.withLongOpt("harmo");
		OPTIONS.addOption(OptionBuilder.create("h"));                
                
	}

	public Configuration(String... args) throws ParseException {
		
		final CommandLine commandLine = new PosixParser().parse(OPTIONS, args, false);
		
		inputPaths = commandLine.getOptionValues('i');
			try {
				if (commandLine.hasOption('I')) {
					inputDataType = RandomAccessGenotypeDataReaderFormats.valueOf(commandLine.getOptionValue('I').toUpperCase());
				} else {
					if (inputPaths[0].endsWith(".vcf")) {
						throw new ParseException("Only vcf.gz is supported. Please see manual on how to do create a vcf.gz file.");
					}
					try {
						inputDataType = RandomAccessGenotypeDataReaderFormats.matchFormatToPath(inputPaths[0]);
					} catch (GenotypeDataException e) {
						throw new ParseException("Unable to determine reference data type based on specified path. Please specify --genotypesType");
					}
				}
			} catch (IllegalArgumentException e) {
				throw new ParseException("Error parsing --genotypesType \"" + commandLine.getOptionValue('I') + "\" is not a valid input data format");
			}
			
		outputFile = new File(commandLine.getOptionValue('o'));
		risksnpsFile = new File(commandLine.getOptionValue('r'));
		onlyCount = commandLine.getOptionValue('c');	
                harmonizedData = commandLine.getOptionValue('h');
		inclusionThreshold = Double.valueOf(commandLine.getOptionValue('t'));		
	}

	public void printOptions() {

		System.out.println("Interpreted arguments: ");


		System.out.print(" - input genotypes " + inputDataType.getName() + ":");
		LOGGER.info("input genotypes " + inputDataType.getName() + ":");
		for (String path : inputPaths) {
			System.out.print(" " + path);
			LOGGER.info(" " + path);
		}
		System.out.println();


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

	public String[] getInputPaths() {
		return inputPaths;
	}

	public RandomAccessGenotypeDataReaderFormats getInputDataType() {
		return inputDataType;
	}

	public File getRisksnpsFile() {
		return risksnpsFile;
	}
	
	public File getOutputFile() {
		return outputFile;
	}	

	public String getOnlyCount() {
		return onlyCount;
	}
	public String getHarmonizedData() {
		return harmonizedData;
	}
        public double getInclusionThreshold() {
		return inclusionThreshold;
	}

}
