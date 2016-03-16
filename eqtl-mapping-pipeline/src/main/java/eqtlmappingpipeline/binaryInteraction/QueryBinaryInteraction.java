package eqtlmappingpipeline.binaryInteraction;

import au.com.bytecode.opencsv.CSVReader;
import au.com.bytecode.opencsv.CSVWriter;
import eqtlmappingpipeline.Main;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.io.Writer;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import umcg.genetica.containers.Pair;
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
		OptionBuilder.withDescription("Gene name (optional)");
		OptionBuilder.withLongOpt("gene");
		OPTIONS.addOption(OptionBuilder.create("g"));

		OptionBuilder.withArgName("string");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Covariate name (optional)");
		OptionBuilder.withLongOpt("covariate");
		OPTIONS.addOption(OptionBuilder.create("c"));

		OptionBuilder.withArgName("string");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Variant name (optional)");
		OptionBuilder.withLongOpt("variant");
		OPTIONS.addOption(OptionBuilder.create("v"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("File with queries. Must have header. All columns are optional, options are gene, variant and covariate. Any combination of headers is possible (optional)");
		OptionBuilder.withLongOpt("queryFile");
		OPTIONS.addOption(OptionBuilder.create("qf"));

		OptionBuilder.withArgName("double");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Minimum absolute interaction z-score (not yet implemented)");
		OptionBuilder.withLongOpt("interactionZ");
		OPTIONS.addOption(OptionBuilder.create("iz"));

		OptionBuilder.withDescription("Only output meta z-scores");
		OptionBuilder.withLongOpt("onlyMetaZ");
		OPTIONS.addOption(OptionBuilder.create("oz"));

	}

	public static void main(String[] args) throws UnsupportedEncodingException, IOException, Exception {

		final File inputInteractionFile;
		final File outputFile;
		final String queryGeneName;
		final String queryCovariateName;
		final String queryVariantName;
		final double queryMinAbsInteractionZ;
		final File queryFile;
		final boolean onlyOutputMetaZ;

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

			if (commandLine.hasOption("qf")) {
				queryFile = new File(commandLine.getOptionValue("qf"));
				if (queryGeneName != null || queryVariantName != null || queryCovariateName != null) {
					System.err.println("Cannot combine query file with commandline query arguments");
					System.exit(1);
					return;
				}
			} else {
				queryFile = null;
			}

			if (commandLine.hasOption("iz")) {
				try {
					queryMinAbsInteractionZ = Double.parseDouble(commandLine.getOptionValue("iz"));
				} catch (NumberFormatException ex) {
					System.err.println("Cannot not parse interactionZ as double: " + commandLine.getOptionValue("iz"));
					System.exit(1);
					return;
				}
			} else {
				queryMinAbsInteractionZ = -1;
			}

			onlyOutputMetaZ = commandLine.hasOption("oz");

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
		if (outputFile != null && !onlyOutputMetaZ) {
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
		if (queryFile != null) {
			outputWriter.write("# - Query file: " + queryFile.getAbsolutePath());
			outputWriter.write('\n');
		}
		if (onlyOutputMetaZ) {
			outputWriter.write("# - Only outputing meta z-scores");
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

		outputWriter.flush();

		final LinkedHashSet<InteractoinQuery> interactionQueries;
		final boolean interactionQueriesOnlyCovariates;
		if (queryFile != null) {
			Pair<LinkedHashSet<InteractoinQuery>, Boolean> loadRes = loadInteractionQueries(queryFile);
			interactionQueries = loadRes.getLeft();
			interactionQueriesOnlyCovariates = loadRes.getRight();
		} else {
			interactionQueries = null;
			interactionQueriesOnlyCovariates = false;
		}

		CSVWriter tableWriter;
		if (onlyOutputMetaZ) {
			if (outputFile == null) {
				throw new Exception("Use of option --onlyMetaZ only possible in combination with output file");
			}
			tableWriter = new CSVWriter(new BufferedWriter(new FileWriter(outputFile)), '\t', '\0', '\0');
		} else {
			tableWriter = new CSVWriter(outputWriter, '\t', '\0', '\0');
		}

		if (onlyOutputMetaZ && !inputFile.isMetaAnalysis()) {
			throw new Exception("No meta analysis information detected cannot use option: --onlyMetaZ");
		}

		final String[] row;

		if (onlyOutputMetaZ) {
			row = new String[1];
		} else {
			final int columnCount = 7
					+ ((5 + (inputFile.isNormalQtlStored() ? 2 : 0) + (inputFile.isFlippedZscoreStored() ? 1 : 0)) * inputFile.getCohortCount())
					+ (inputFile.isMetaAnalysis() ? (3 + (inputFile.isNormalQtlStored() ? 1 : 0) + (inputFile.isFlippedZscoreStored() ? 1 : 0)) : 0);



			row = new String[columnCount];
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
		}

		if (interactionQueries != null) {
			
			if (interactionQueriesOnlyCovariates) {
				HashSet<String> covariateNames = new HashSet<>(interactionQueries.size());
				for (InteractoinQuery interactionQuery : interactionQueries) {
					covariateNames.add(interactionQuery.getCovariate());
				}
				doQueryCovariates(covariateNames, inputFile, tableWriter, row, onlyOutputMetaZ);
			} else {
				for (InteractoinQuery interactionQuery : interactionQueries) {
					doQuery(interactionQuery.getGene(), interactionQuery.getVariant(), interactionQuery.getCovariate(), inputFile, tableWriter, row, onlyOutputMetaZ);
				}
			}
			
		} else {
			doQuery(queryGeneName, queryVariantName, queryCovariateName, inputFile, tableWriter, row, onlyOutputMetaZ);
		}



		tableWriter.close();
		outputWriter.close();

	}

	@SuppressWarnings({"null", "ConstantConditions"})
	private static void addRow(BinaryInteractionQueryResult queryRestult, BinaryInteractionFile inputFile, CSVWriter tableWriter, String[] row, boolean onlyOutputMetaZ) throws BinaryInteractionFileException, IOException {

		if (onlyOutputMetaZ) {
			BinaryInteractionZscores zscroresInteraction = queryRestult.getInteractionZscores();
			row[0] = String.valueOf(zscroresInteraction.getZscoreInteractionMeta());
		} else {
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
		}


		tableWriter.writeNext(row);
	}

	private static void doQuery(final String queryGeneName, final String queryVariantName, final String queryCovariateName, BinaryInteractionFile inputFile, CSVWriter tableWriter, String[] row, boolean onlyOutputMetaZ) throws IOException, BinaryInteractionFileException {

		if (queryGeneName != null && queryVariantName != null && queryCovariateName != null) {

			addRow(inputFile.readVariantGeneCovariateResults(queryVariantName, queryGeneName, queryCovariateName), inputFile, tableWriter, row, onlyOutputMetaZ);

		} else if (queryGeneName != null && queryVariantName != null) {

			for (Iterator<BinaryInteractionQueryResult> iterator = inputFile.readVariantGeneResults(queryVariantName, queryGeneName); iterator.hasNext();) {
				addRow(iterator.next(), inputFile, tableWriter, row, onlyOutputMetaZ);
			}

		} else if (queryVariantName != null) {

			int[] genePointers = inputFile.getVariant(queryVariantName).getGenePointers();
			for (int genePointer : genePointers) {

				BinaryInteractionGene gene = inputFile.getGene(genePointer);
				if (queryCovariateName != null) {

					if (inputFile.containsInteraction(queryVariantName, gene.getName(), queryCovariateName)) {
						addRow(inputFile.readVariantGeneCovariateResults(queryVariantName, gene.getName(), queryCovariateName), inputFile, tableWriter, row, onlyOutputMetaZ);
					}

				} else {
					for (Iterator<BinaryInteractionQueryResult> iterator = inputFile.readVariantGeneResults(queryVariantName, gene.getName()); iterator.hasNext();) {
						addRow(iterator.next(), inputFile, tableWriter, row, onlyOutputMetaZ);
					}
				}

			}

		} else if (queryGeneName != null) {

			int[] variantPointers = inputFile.getGene(queryGeneName).getVariantPointers();
			for (int variantPointer : variantPointers) {

				BinaryInteractionVariant variant = inputFile.getVariant(variantPointer);
				if (queryCovariateName != null) {

					if (inputFile.containsInteraction(variant.getName(), queryGeneName, queryCovariateName)) {
						addRow(inputFile.readVariantGeneCovariateResults(variant.getName(), queryGeneName, queryCovariateName), inputFile, tableWriter, row, onlyOutputMetaZ);
					}

				} else {
					for (Iterator<BinaryInteractionQueryResult> iterator = inputFile.readVariantGeneResults(variant.getName(), queryGeneName); iterator.hasNext();) {
						addRow(iterator.next(), inputFile, tableWriter, row, onlyOutputMetaZ);
					}
				}


			}

		} else {

			for (BinaryInteractionVariant variant : inputFile.getVariants()) {

				String variantName = variant.getName();

				int[] genePointers = inputFile.getVariant(variantName).getGenePointers();
				for (int genePointer : genePointers) {

					BinaryInteractionGene gene = inputFile.getGene(genePointer);
					if (queryCovariateName != null) {

						if (inputFile.containsInteraction(variantName, gene.getName(), queryCovariateName)) {
							addRow(inputFile.readVariantGeneCovariateResults(variantName, gene.getName(), queryCovariateName), inputFile, tableWriter, row, onlyOutputMetaZ);
						}

					} else {
						for (Iterator<BinaryInteractionQueryResult> iterator = inputFile.readVariantGeneResults(variantName, gene.getName()); iterator.hasNext();) {
							addRow(iterator.next(), inputFile, tableWriter, row, onlyOutputMetaZ);
						}
					}

				}

			}



		}
	}

	private static void doQueryCovariates(final HashSet<String> queryCovariateNames, BinaryInteractionFile inputFile, CSVWriter tableWriter, String[] row, boolean onlyOutputMetaZ) throws IOException, BinaryInteractionFileException {

		for (BinaryInteractionVariant variant : inputFile.getVariants()) {

			String variantName = variant.getName();

			int[] genePointers = inputFile.getVariant(variantName).getGenePointers();
			for (int genePointer : genePointers) {

				BinaryInteractionGene gene = inputFile.getGene(genePointer);

				for (Iterator<BinaryInteractionQueryResult> iterator = inputFile.readVariantGeneResults(variantName, gene.getName()); iterator.hasNext();) {
					BinaryInteractionQueryResult next = iterator.next();
					if(queryCovariateNames.contains(next.getCovariateName())){
						addRow(next, inputFile, tableWriter, row, onlyOutputMetaZ);
					}
					
				}

			}

		}

	}

	private static Pair<LinkedHashSet<InteractoinQuery>, Boolean> loadInteractionQueries(File queryFile) throws FileNotFoundException, IOException, Exception {

		LinkedHashSet<InteractoinQuery> interactionQueries = new LinkedHashSet<InteractoinQuery>();
		final CSVReader queryReader = new CSVReader(new FileReader(queryFile), '\t', '\0');

		String[] nextLine = queryReader.readNext();

		int variantCol = -1;
		int geneCol = -1;
		int covariateCol = -1;

		//Parse header
		for (int i = 0; i < nextLine.length; ++i) {
			String headerEntry = nextLine[i].toLowerCase();
			switch (headerEntry) {
				case "variant":
					if (variantCol != -1) {
						throw new Exception("Variant column found twice");
					}
					variantCol = i;
					break;
				case "gene":
					if (geneCol != -1) {
						throw new Exception("Gene column found twice");
					}
					geneCol = i;
					break;
				case "covariate":
					if (covariateCol != -1) {
						throw new Exception("Covariate column found twice");
					}
					covariateCol = i;
					break;

			}

		}

		if (variantCol == -1 && geneCol == -1 && covariateCol == -1) {
			throw new Exception("Did not detect appropiate header in query file");

		}

		while ((nextLine = queryReader.readNext()) != null) {
			String variant = null;
			String gene = null;
			String covariate = null;

			if (variantCol != -1) {
				variant = nextLine[variantCol];
			}
			if (geneCol != -1) {
				gene = nextLine[geneCol];
			}
			if (covariateCol != -1) {
				covariate = nextLine[covariateCol];
			}
			interactionQueries.add(new InteractoinQuery(variant, gene, covariate));
		}
		queryReader.close();

		return new Pair(interactionQueries, variantCol == -1 && geneCol == -1);
	}

	private static class InteractoinQuery {

		private final String variant;
		private final String gene;
		private final String covariate;

		public InteractoinQuery(String variant, String gene, String covariate) {
			this.variant = variant;
			this.gene = gene;
			this.covariate = covariate;
		}

		public String getVariant() {
			return variant;
		}

		public String getGene() {
			return gene;
		}

		public String getCovariate() {
			return covariate;
		}

		@Override
		public int hashCode() {
			int hash = 5;
			hash = 67 * hash + (this.variant != null ? this.variant.hashCode() : 0);
			hash = 67 * hash + (this.gene != null ? this.gene.hashCode() : 0);
			hash = 67 * hash + (this.covariate != null ? this.covariate.hashCode() : 0);
			return hash;
		}

		@Override
		public boolean equals(Object obj) {
			if (obj == null) {
				return false;
			}
			if (getClass() != obj.getClass()) {
				return false;
			}
			final InteractoinQuery other = (InteractoinQuery) obj;
			if ((this.variant == null) ? (other.variant != null) : !this.variant.equals(other.variant)) {
				return false;
			}
			if ((this.gene == null) ? (other.gene != null) : !this.gene.equals(other.gene)) {
				return false;
			}
			if ((this.covariate == null) ? (other.covariate != null) : !this.covariate.equals(other.covariate)) {
				return false;
			}
			return true;
		}
	}
}
