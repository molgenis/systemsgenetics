package eqtlmappingpipeline.binaryInteraction;

import au.com.bytecode.opencsv.CSVWriter;
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
import java.util.Iterator;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import umcg.genetica.io.binInteraction.BinaryInteractionCohort;
import umcg.genetica.io.binInteraction.BinaryInteractionFile;
import umcg.genetica.io.binInteraction.BinaryInteractionFileException;
import umcg.genetica.io.binInteraction.BinaryInteractionQtlZscores;
import umcg.genetica.io.binInteraction.BinaryInteractionQueryResult;
import umcg.genetica.io.binInteraction.BinaryInteractionZscores;
import umcg.genetica.io.binInteraction.gene.BinaryInteractionGene;
import umcg.genetica.io.binInteraction.variant.BinaryInteractionVariant;

/**
 *
 * @author Patrick Deelen
 */
public class QueryBinaryInteraction {

	private static final String VERSION = Main.VERSION;
	private static final DateFormat DATE_TIME_FORMAT = new SimpleDateFormat("dd-MM-yyyy HH:mm:ss");
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
		OptionBuilder.withDescription("Output file (optional)");
		OptionBuilder.withLongOpt("output");
		OPTIONS.addOption(OptionBuilder.create("o"));

		OptionBuilder.withArgName("string");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Gene name");
		OptionBuilder.withLongOpt("gene");
		OPTIONS.addOption(OptionBuilder.create("g"));

		OptionBuilder.withArgName("string");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Covariate name (optional)");
		OptionBuilder.withLongOpt("cocariate");
		OPTIONS.addOption(OptionBuilder.create("c"));

		OptionBuilder.withArgName("string");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Variant name");
		OptionBuilder.withLongOpt("variant");
		OPTIONS.addOption(OptionBuilder.create("v"));

		OptionBuilder.withArgName("double");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Minimum absolute interaction z-score (not yet implemented)");
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
			if (commandLine.hasOption("o")) {
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
		if (outputFile != null) {
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
		if (outputFile != null) {
			outputWriter.write("# - Output file: " + outputFile.getAbsolutePath());
			outputWriter.write('\n');
		}
		if (queryGeneName != null) {
			outputWriter.write("# - Query gene: " + queryGeneName);
			outputWriter.write('\n');
		}
		if (queryCovariateName != null) {
			outputWriter.write("# - Query covariate: " + queryCovariateName);
			outputWriter.write('\n');
		}
		if (queryVariantName != null) {
			outputWriter.write("# - Query variant: " + queryVariantName);
			outputWriter.write('\n');
		}
		if (queryMinAbsInteractionZ > 0) {
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
		for (BinaryInteractionCohort cohort : inputFile.getCohorts()) {
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



		CSVWriter tableWriter = new CSVWriter(outputWriter, '\t', '\0', '\0');

		int columnCount =
				7
				+ ((5 + (inputFile.isNormalQtlStored() ? 2 : 0) + (inputFile.isFlippedZscoreStored() ? 1 : 0)) * inputFile.getCohortCount())
				+ (inputFile.isMetaAnalysis() ? (3 + (inputFile.isNormalQtlStored() ? 1 : 0) + (inputFile.isFlippedZscoreStored() ? 1 : 0)) : 0);


		String[] row = new String[columnCount];
		int c = 0;

		row[c++] = "Variant";
		row[c++] = "Gene";
		row[c++] = "Covariate";
		row[c++] = "Variant_chr";
		row[c++] = "Variant_pos";
		row[c++] = "Variant alleles";
		row[c++] = "Assessed_allele";

		for (BinaryInteractionCohort cohort : inputFile.getCohorts()) {

			String cohortName = cohort.getName();

			if (inputFile.isNormalQtlStored()) {
				row[c++] = cohortName + "_QTL_sample_count";
				row[c++] = cohortName + "_QTL_Z-score";
			}

			row[c++] = cohortName + "_interaction_sample_count";
			row[c++] = cohortName + "_interaction_r2";
			row[c++] = cohortName + "_variant_Z-score";
			row[c++] = cohortName + "_covariate_Z-score";
			row[c++] = cohortName + "_interaction_Z-score";

			if (inputFile.isFlippedZscoreStored()) {
				row[c++] = cohortName + "_flipped_interaction_Z-score";
			}

		}

		if (inputFile.isMetaAnalysis()) {
			if (inputFile.isNormalQtlStored()) {
				row[c++] = "Meta_QTL_Z-score";
			}
			row[c++] = "Meta_variant_Z-score";
			row[c++] = "Meta_covariate_Z-score";
			row[c++] = "Meta_interaction_Z-score";
			if (inputFile.isFlippedZscoreStored()) {
				row[c++] = "Meta_flipped_interaction_Z-score";
			}
		}

		tableWriter.writeNext(row);

		if (queryGeneName != null && queryVariantName != null && queryCovariateName != null) {

			addRow(inputFile.readVariantGeneCovariateResults(queryVariantName, queryGeneName, queryCovariateName), inputFile, tableWriter, row);

		} else if (queryGeneName != null && queryVariantName != null) {

			for (Iterator<BinaryInteractionQueryResult> iterator = inputFile.readVariantGeneResults(queryVariantName, queryGeneName); iterator.hasNext();) {
				addRow(iterator.next(), inputFile, tableWriter, row);
			}

		} else if (queryVariantName != null) {

			int[] genePointers = inputFile.getVariant(queryVariantName).getGenePointers();
			for (int genePointer : genePointers) {

				BinaryInteractionGene gene = inputFile.getGene(genePointer);
				if (queryCovariateName != null) {

					addRow(inputFile.readVariantGeneCovariateResults(queryVariantName, gene.getName(), queryCovariateName), inputFile, tableWriter, row);

				} else {
					for (Iterator<BinaryInteractionQueryResult> iterator = inputFile.readVariantGeneResults(queryVariantName, gene.getName()); iterator.hasNext();) {
						addRow(iterator.next(), inputFile, tableWriter, row);
					}
				}


			}

		} else {
			tableWriter.close();
			outputWriter.append("ERROR not yet supported");
			System.err.println("ERROR not yet supported");
		}

		tableWriter.close();
		outputWriter.close();

	}

	@SuppressWarnings({"null", "ConstantConditions"})
	private static void addRow(BinaryInteractionQueryResult queryRestult, BinaryInteractionFile inputFile, CSVWriter tableWriter, String[] row) throws BinaryInteractionFileException, IOException {
		int c = 0;

		row[c++] = queryRestult.getVariantName();
		row[c++] = queryRestult.getGeneName();
		row[c++] = queryRestult.getCovariateName();

		BinaryInteractionVariant variant = inputFile.getVariant(queryRestult.getVariantName());
		row[c++] = variant.getChr();
		row[c++] = String.valueOf(variant.getPos());
		row[c++] = variant.getRefAllele().getAlleleAsString() + '/' + variant.getAltAllele().getAlleleAsString();
		row[c++] = variant.getAltAllele().toString();

		BinaryInteractionQtlZscores zscroresQtl = queryRestult.getQtlZscores();
		BinaryInteractionZscores zscroresInteraction = queryRestult.getInteractionZscores();

		for (int cohortIndex = 0; cohortIndex < inputFile.getCohortCount(); ++cohortIndex) {

			if (inputFile.isNormalQtlStored()) {
				row[c++] = String.valueOf(zscroresQtl.getSampleCounts()[cohortIndex]);
				row[c++] = String.valueOf(zscroresQtl.getZscores()[cohortIndex]);
			}

			row[c++] = String.valueOf(zscroresInteraction.getSamplesInteractionCohort()[cohortIndex]);
			row[c++] = String.valueOf(zscroresInteraction.getrSquaredCohort()[cohortIndex]);
			row[c++] = String.valueOf(zscroresInteraction.getZscoreSnpCohort()[cohortIndex]);
			row[c++] = String.valueOf(zscroresInteraction.getZscoreCovariateCohort()[cohortIndex]);
			row[c++] = String.valueOf(zscroresInteraction.getZscoreInteractionCohort()[cohortIndex]);

			if (inputFile.isFlippedZscoreStored()) {
				row[c++] = String.valueOf(zscroresInteraction.getZscoreInteractionFlippedCohort()[cohortIndex]);
			}

		}

		if (inputFile.isMetaAnalysis()) {
			if (inputFile.isNormalQtlStored()) {
				row[c++] = String.valueOf(zscroresQtl.getMetaZscore());
			}
			row[c++] = String.valueOf(zscroresInteraction.getZscoreSnpMeta());
			row[c++] = String.valueOf(zscroresInteraction.getZscoreCovariateMeta());
			row[c++] = String.valueOf(zscroresInteraction.getZscoreInteractionMeta());
			if (inputFile.isFlippedZscoreStored()) {
				row[c++] = String.valueOf(zscroresInteraction.getZscoreInteractionFlippedMeta());
			}
		}

		tableWriter.writeNext(row);
	}
}
