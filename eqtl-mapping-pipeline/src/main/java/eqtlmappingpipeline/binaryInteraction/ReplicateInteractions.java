package eqtlmappingpipeline.binaryInteraction;

import au.com.bytecode.opencsv.CSVWriter;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DateFormat;
import java.text.NumberFormat;
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
public class ReplicateInteractions {

	private static final DateFormat DATE_TIME_FORMAT = new SimpleDateFormat("dd-MM-yyyy HH:mm:ss");
	private static final Date currentDataTime = new Date();
	private static final Options OPTIONS;

	static {

		OPTIONS = new Options();

		Option option;

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Binary interaction file (must be a meta analysis)");
		OptionBuilder.withLongOpt("input");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("i"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Binary interaction file to use as replication (must be a meta analysis)");
		OptionBuilder.withLongOpt("replication");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("r"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Ouput prefix");
		OptionBuilder.withLongOpt("output");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("o"));

		OptionBuilder.withArgName("double");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Minimum absolute interaction z-score");
		OptionBuilder.withLongOpt("interactionZ");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("iz"));

		OptionBuilder.withArgName("double");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Minimum absolute replication interaction z-score");
		OptionBuilder.withLongOpt("replicationInteractionZ");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("riz"));

		OptionBuilder.withDescription("If set match variant on chr-pos");
		OptionBuilder.withLongOpt("chrPos");
		OPTIONS.addOption(OptionBuilder.create("cp"));

	}

	public static void main(String[] args) throws FileNotFoundException, IOException, BinaryInteractionFileException {

		final File inputInteractionFile;
		final File replicationInteractionFile;
		final double minAbsInteractionZ;
		final double minAbsReplicationInteractionZ;
		final boolean matchOnChrPos;
		final String outputPrefix;

		try {
			final CommandLine commandLine = new PosixParser().parse(OPTIONS, args, false);

			inputInteractionFile = new File(commandLine.getOptionValue("i"));
			replicationInteractionFile = new File(commandLine.getOptionValue("r"));
			outputPrefix = commandLine.getOptionValue("o");

			try {
				minAbsInteractionZ = Double.parseDouble(commandLine.getOptionValue("iz"));
			} catch (NumberFormatException ex) {
				System.out.println("Cannot not parse --interactionZ as double: " + commandLine.getOptionValue("iz"));
				System.exit(1);
				return;
			}

			try {
				minAbsReplicationInteractionZ = Double.parseDouble(commandLine.getOptionValue("riz"));
			} catch (NumberFormatException ex) {
				System.out.println("Cannot not parse --replicationInteractionZ as double: " + commandLine.getOptionValue("riz"));
				System.exit(1);
				return;
			}

			matchOnChrPos = commandLine.hasOption("cp");

		} catch (ParseException ex) {
			System.err.println("Invalid command line arguments: ");
			System.err.println(ex.getMessage());
			System.err.println();
			new HelpFormatter().printHelp(" ", OPTIONS);
			System.exit(1);
			return;
		}

		System.out.println("Input file: " + inputInteractionFile.getAbsolutePath());
		System.out.println("Replication file: " + replicationInteractionFile.getAbsolutePath());
		System.out.println("Output prefix: " + outputPrefix);
		System.out.println("Min interaction z-score: " + minAbsInteractionZ);
		System.out.println("Min replication interaction z-score: " + minAbsReplicationInteractionZ);
		if (matchOnChrPos) {
			System.out.println("Matching variants on chr-pos");
		}
		System.out.println("");

		BinaryInteractionFile inputFile = BinaryInteractionFile.load(inputInteractionFile, true);
		BinaryInteractionFile replicationFile = BinaryInteractionFile.load(replicationInteractionFile, true);

		CSVWriter replicatedSameDirectionWriter = new CSVWriter(new BufferedWriter(new FileWriter(outputPrefix + "replicatedSameDirection.txt")), '\t', '\0', '\0');

		String[] row = new String[15];
		int c = 0;

		row[c++] = "Variant";
		row[c++] = "Gene";
		row[c++] = "Covariate";
		row[c++] = "Variant_chr";
		row[c++] = "Variant_pos";
		row[c++] = "Variant alleles";
		row[c++] = "Assessed_allele";
		row[c++] = "Discovery_meta_QTL";
		row[c++] = "Discovery_meta_SNP";
		row[c++] = "Discovery_meta_covariate";
		row[c++] = "Discovery_meta_interaction";
		row[c++] = "Replication_meta_QTL";
		row[c++] = "Replication_meta_SNP";
		row[c++] = "Replication_meta_covariate";
		row[c++] = "Replication_meta_interaction";

		replicatedSameDirectionWriter.writeNext(row);


		int significant = 0;
		int notSignificant = 0;
		int notTestedInReplication = 0;
		int notSignificantReplicationSameDirection = 0;
		int notSignificantReplicationOppositeDirection = 0;
		int significantReplicationOppositeDirection = 0;
		int significantReplicationSameDirection = 0;

		int reporter = 0;

		for (BinaryInteractionVariant variant : inputFile.getVariants()) {

			String variantName = variant.getName();

			BinaryInteractionVariant replicationVariant;

			if (matchOnChrPos) {
				replicationVariant = replicationFile.getVariant(variant.getChr(), variant.getPos());
			} else {
				if (replicationFile.containsVariant(variantName)) {
					replicationVariant = replicationFile.getVariant(variantName);
				} else {
					replicationVariant = null;
				}
			}
			//Do loop anyway to also count not replicated

			int[] genePointers = inputFile.getVariant(variantName).getGenePointers();
			for (int genePointer : genePointers) {

				BinaryInteractionGene gene = inputFile.getGene(genePointer);

				covairates:
				for (Iterator<BinaryInteractionQueryResult> iterator = inputFile.readVariantGeneResults(variantName, gene.getName()); iterator.hasNext();) {

					BinaryInteractionQueryResult interation = iterator.next();

					double metaInteractionZ = interation.getInteractionZscores().getZscoreInteractionMeta();

					if (metaInteractionZ >= minAbsInteractionZ || metaInteractionZ <= -minAbsInteractionZ) {
						++significant;

						if (replicationVariant != null && replicationFile.containsInteraction(replicationVariant.getName(), gene.getName(), interation.getCovariateName())) {

							if (!(variant.getRefAllele() == replicationVariant.getRefAllele() && variant.getAltAllele() == replicationVariant.getAltAllele())
									&& !(variant.getRefAllele() == replicationVariant.getAltAllele() && variant.getAltAllele() == replicationVariant.getRefAllele())) {
								System.err.println("Allele mismatch!");
								continue covairates;
							}

							BinaryInteractionZscores replicationZscores = replicationFile.readInteractionResults(replicationVariant.getName(), gene.getName(), interation.getCovariateName());
							double replicationInteractionZscore = replicationZscores.getZscoreInteractionMeta();

							if (variant.getAltAllele() != replicationVariant.getAltAllele()) {
								replicationInteractionZscore *= -1;
							}

							if (replicationInteractionZscore <= -minAbsReplicationInteractionZ || replicationInteractionZscore >= minAbsReplicationInteractionZ) {
								if (metaInteractionZ * replicationInteractionZscore >= 0) {
									++significantReplicationSameDirection;
									
									BinaryInteractionQtlZscores replicationQtlRes = replicationFile.readQtlResults(replicationVariant.getName(), gene.getName());

									c = 0;
									row[c++] = variantName;
									row[c++] = gene.getName();
									row[c++] = interation.getCovariateName();
									row[c++] = variant.getChr();
									row[c++] = String.valueOf(variant.getPos());
									row[c++] = variant.getRefAllele().getAlleleAsString() + "/" + variant.getAltAllele().getAlleleAsString();
									row[c++] = variant.getAltAllele().getAlleleAsString();
									row[c++] = String.valueOf(interation.getQtlZscores().getMetaZscore());
									row[c++] = String.valueOf(interation.getInteractionZscores().getZscoreSnpMeta());
									row[c++] = String.valueOf(interation.getInteractionZscores().getZscoreCovariateMeta());
									row[c++] = String.valueOf(interation.getInteractionZscores().getZscoreInteractionMeta());
									row[c++] = String.valueOf(replicationQtlRes.getMetaZscore());
									row[c++] = String.valueOf(replicationZscores.getZscoreSnpMeta());
									row[c++] = String.valueOf(replicationZscores.getZscoreCovariateMeta());
									row[c++] = String.valueOf(replicationZscores.getZscoreInteractionMeta());

									replicatedSameDirectionWriter.writeNext(row);

								} else {
									++significantReplicationOppositeDirection;
								}
							} else {
								if (metaInteractionZ * replicationInteractionZscore >= 0) {
									++notSignificantReplicationSameDirection;
								} else {
									++notSignificantReplicationOppositeDirection;
								}
							}


						} else {
							++notTestedInReplication;
						}

					} else {
						++notSignificant;
					}
				}

				++reporter;
				if (reporter % 500 == 0) {
					System.out.println("Parsed " + reporter + " of " + inputFile.getVariantGeneCombinations() + " variant-gene combinations");
				}

			}

		}

		replicatedSameDirectionWriter.close();

		NumberFormat numberFormat = NumberFormat.getInstance();
		numberFormat.setMinimumFractionDigits(0);
		numberFormat.setMaximumFractionDigits(2);

		System.out.println("");
		System.out.println("Total number of interactions: " + numberFormat.format(notSignificant + significant));
		System.out.println(" - Not significant: " + numberFormat.format(notSignificant) + " (" + numberFormat.format(notSignificant * 100d / (notSignificant + significant)) + "%)");
		System.out.println(" - Significant: " + numberFormat.format(significant) + " (" + numberFormat.format(significant * 100d / (notSignificant + significant)) + "%)");
		System.out.println("  * Not in replication: " + numberFormat.format(notTestedInReplication) + " (" + numberFormat.format(notTestedInReplication * 100d / significant) + "%)");
		System.out.println("  * Not significant in replication: " + numberFormat.format(notSignificantReplicationSameDirection + notSignificantReplicationOppositeDirection) + " (" + numberFormat.format((notSignificantReplicationSameDirection + notSignificantReplicationOppositeDirection) * 100d / significant) + "%)");
		System.out.println("   # Same direction: " + numberFormat.format(notSignificantReplicationSameDirection) + " (" + numberFormat.format(notSignificantReplicationSameDirection * 100d / (notSignificantReplicationSameDirection + notSignificantReplicationOppositeDirection)) + "%)");
		System.out.println("   # Opposite direction: " + numberFormat.format(notSignificantReplicationOppositeDirection) + " (" + numberFormat.format(notSignificantReplicationOppositeDirection * 100d / (notSignificantReplicationSameDirection + notSignificantReplicationOppositeDirection)) + "%)");
		System.out.println("  * Significant in replication: " + numberFormat.format(significantReplicationSameDirection + significantReplicationOppositeDirection) + " (" + numberFormat.format((significantReplicationSameDirection + significantReplicationOppositeDirection) * 100d / significant) + "%)");
		System.out.println("   # Same direction: " + numberFormat.format(significantReplicationSameDirection) + " (" + numberFormat.format(significantReplicationSameDirection * 100d / (significantReplicationSameDirection + significantReplicationOppositeDirection)) + "%)");
		System.out.println("   # Opposite direction: " + numberFormat.format(significantReplicationOppositeDirection) + " (" + numberFormat.format(significantReplicationOppositeDirection * 100d / (significantReplicationSameDirection + significantReplicationOppositeDirection)) + "%)");

	}
}
