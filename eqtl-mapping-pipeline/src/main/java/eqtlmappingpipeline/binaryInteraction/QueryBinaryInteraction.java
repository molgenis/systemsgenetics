package eqtlmappingpipeline.binaryInteraction;

import eqtlmappingpipeline.Main;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.io.Writer;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import umcg.genetica.io.binInteraction.BinaryInteractionCohort;
import umcg.genetica.io.binInteraction.BinaryInteractionFile;
import umcg.genetica.io.binInteraction.BinaryInteractionQtlZscores;
import umcg.genetica.io.binInteraction.BinaryInteractionZscores;
import umcg.genetica.io.binInteraction.variant.BinaryInteractionVariant;

/**
 *
 * @author Patrick Deelen
 */
public class QueryBinaryInteraction {

	private static final String VERSION = Main.VERSION;
	private static final DateFormat DATE_TIME_FORMAT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
	private static final Date currentDataTime = new Date();
	private static final Options OPTIONS;

	static {

		OPTIONS = new Options();

		Option option;

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Binary interaction file");
		OptionBuilder.withLongOpt("input");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("i"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Output file");
		OptionBuilder.withLongOpt("output");
		OPTIONS.addOption(OptionBuilder.create("o"));

		OptionBuilder.withArgName("string");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Gene name");
		OptionBuilder.withLongOpt("gene");
		OPTIONS.addOption(OptionBuilder.create("g"));

		OptionBuilder.withArgName("string");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Covariate name");
		OptionBuilder.withLongOpt("cocariate");
		OPTIONS.addOption(OptionBuilder.create("c"));

		OptionBuilder.withArgName("string");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Variant name");
		OptionBuilder.withLongOpt("variant");
		OPTIONS.addOption(OptionBuilder.create("v"));

		OptionBuilder.withArgName("string");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Minimum absolute interaction z-score");
		OptionBuilder.withLongOpt("interactionZ");
		OPTIONS.addOption(OptionBuilder.create("iz"));
		
	}

	public static void main(String[] args) throws UnsupportedEncodingException, IOException, Exception {

		final File inputInteractionFile;
		final File outputFile;
		final String queryGeneName;
		final String queryCovariateName;
		final String queryVariantName;
		final double queryMinAbsInteractionZ;

		try {
			final CommandLine commandLine = new PosixParser().parse(OPTIONS, args, false);

			inputInteractionFile = new File(commandLine.getOptionValue("i"));
			if(commandLine.hasOption("o")){
				outputFile = new File(commandLine.getOptionValue("o"));
			} else {
				outputFile = null;
			}
			queryGeneName = commandLine.getOptionValue("g");
			queryCovariateName = commandLine.getOptionValue("c");
			queryVariantName = commandLine.getOptionValue("v");

			if (commandLine.hasOption("iz")) {
				try {
					queryMinAbsInteractionZ = Double.parseDouble(commandLine.getOptionValue("iz"));
				} catch (NumberFormatException ex) {
					System.out.println("Cannot not parse interactionZ as double: " + commandLine.getOptionValue("iz"));
					System.exit(1);
					return;
				}
			} else {
				queryMinAbsInteractionZ = -1;
			}

		} catch (ParseException ex) {
			System.err.println("Invalid command line arguments: ");
			System.err.println(ex.getMessage());
			System.err.println();
			new HelpFormatter().printHelp(" ", OPTIONS);
			System.exit(1);
			return;
		}
		
		BinaryInteractionFile inputFile = BinaryInteractionFile.load(inputInteractionFile, true);
		
		final Writer outputWriter;
		if(outputFile != null){
			outputWriter = new BufferedWriter(new FileWriter(outputFile));
		} else {
			outputWriter = new OutputStreamWriter(System.out);
		}
		
		outputWriter.write("# Query result binary interaction file using software version: " + VERSION);
		outputWriter.write('\n');
		outputWriter.write("# Current data and time: " + DATE_TIME_FORMAT.format(currentDataTime));
		outputWriter.write('\n');
		outputWriter.write("# Command line options: ");
		outputWriter.write('\n');
		outputWriter.write("# - Input file: " + inputInteractionFile.getAbsolutePath());
		outputWriter.write('\n');
		if(outputFile != null){
			outputWriter.write("# - Output file: " + outputFile.getAbsolutePath());
			outputWriter.write('\n');
		}
		if(queryGeneName != null){
			outputWriter.write("# - Query gene: " + queryGeneName);
			outputWriter.write('\n');	
		}
		if(queryCovariateName != null){
			outputWriter.write("# - Query covariate: " + queryCovariateName);
			outputWriter.write('\n');	
		}
		if(queryVariantName != null){
			outputWriter.write("# - Query variant: " + queryVariantName);
			outputWriter.write('\n');	
		}
		if(queryMinAbsInteractionZ > 0){
			outputWriter.write("# - Query minimum absote interaction z-score: " + queryMinAbsInteractionZ);
			outputWriter.write('\n');	
		}
		outputWriter.write("#\n");
		
		outputWriter.write("# Interaction file meta data: ");
		outputWriter.write('\n');
		outputWriter.write("# - Description: " + inputFile.getFileDescription());
		outputWriter.write('\n');
		outputWriter.write("# - Creation data: " + inputFile.getCreationDataTimeString());
		outputWriter.write('\n');
		outputWriter.write("# - Cohorts: " + inputFile.getCohortCount());
		outputWriter.write('\n');
		for(BinaryInteractionCohort cohort : inputFile.getCohorts()){
			outputWriter.write("#    * " + cohort.getName() + " (" + cohort.getSampleCount() + ")");
			outputWriter.write('\n');
		}
		outputWriter.write("# - Variants: " + inputFile.getVariantCount());
		outputWriter.write('\n');
		outputWriter.write("# - Genes: " + inputFile.getGeneCount());
		outputWriter.write('\n');
		outputWriter.write("# - Covariats: " + inputFile.getCovariateCount());
		outputWriter.write('\n');
		outputWriter.write("#\n");
		
		
		
		//CSVWriter writer = new CSVWriter(outputWriter, '\t', '\0', '\0');
		
		if(queryGeneName != null && queryVariantName != null && queryCovariateName != null){
						
			BinaryInteractionVariant variant = inputFile.getVariant(queryVariantName);
			outputWriter.write("Chr variant: " + variant.getChr());
			outputWriter.write('\n');
			outputWriter.write("Pos variant: " + variant.getPos());
			outputWriter.write('\n');
			
			BinaryInteractionQtlZscores zscroresQtl = inputFile.readQtlResults(queryVariantName, queryGeneName);
			outputWriter.write("QTL: " + zscroresQtl.getZscores()[0]);
			outputWriter.write('\n');
						
			BinaryInteractionZscores zscroresInteraction = inputFile.readInteractionResults(queryVariantName, queryGeneName, queryCovariateName);
			
			outputWriter.write("SNP: " + zscroresInteraction.getZscoreSnpCohort()[0]);
			outputWriter.write('\n');
			outputWriter.write("Covariate: " + zscroresInteraction.getZscoreCovariateCohort()[0]);
			outputWriter.write('\n');
			outputWriter.write("Interaction: " + zscroresInteraction.getZscoreInteractionCohort()[0]);
			outputWriter.write('\n');
			outputWriter.write("Samples: " + zscroresInteraction.getSamplesInteractionCohort()[0]);
			outputWriter.write('\n');
		}
	
		outputWriter.close();

	}
}
