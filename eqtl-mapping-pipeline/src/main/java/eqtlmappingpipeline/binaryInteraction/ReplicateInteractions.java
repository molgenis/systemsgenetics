package eqtlmappingpipeline.binaryInteraction;

import au.com.bytecode.opencsv.CSVWriter;
import eqtlmappingpipeline.Main;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Writer;
import java.text.NumberFormat;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
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

	private static final Options OPTIONS;

	static {

		OPTIONS = new Options();

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

		OptionBuilder.withArgName("double");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Minimum absolute interaction z-score to count covariate");
		OptionBuilder.withLongOpt("covariateInteractionZ");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("ciz"));

		OptionBuilder.withArgName("double");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Minimum absolute replication interaction z-score to count covariate");
		OptionBuilder.withLongOpt("covariateReplicationInteractionZ");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("criz"));

		OptionBuilder.withDescription("If set match variant on chr-pos");
		OptionBuilder.withLongOpt("chrPos");
		OPTIONS.addOption(OptionBuilder.create("cp"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("File with covariates to include in analysis");
		OptionBuilder.withLongOpt("covariats");
		OPTIONS.addOption(OptionBuilder.create("c"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("File with eQTL genes to include in analysis");
		OptionBuilder.withLongOpt("genes");
		OPTIONS.addOption(OptionBuilder.create("g"));
		
		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("File with covariates to test for replication");
		OptionBuilder.withLongOpt("covariatsReplication");
		OPTIONS.addOption(OptionBuilder.create("cr"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("File with eQTL genes to test for replication");
		OptionBuilder.withLongOpt("genesReplication");
		OPTIONS.addOption(OptionBuilder.create("gr"));

	}

	public static void main(String[] args) throws FileNotFoundException, IOException, BinaryInteractionFileException {

		final File inputInteractionFile;
		final File replicationInteractionFile;
		final double minAbsInteractionZ;
		final double minAbsReplicationInteractionZ;
		final double minAbsInteractionZCovariateCount;
		final double minAbsReplicationInteractionZCovariateCount;
		final boolean matchOnChrPos;
		final String outputPrefix;
		final File covariatesToIncludeFile;
		final File genesToIncludeFile;
		final File covariatesReplicationToIncludeFile;
		final File genesReplicationToIncludeFile;

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

			try {
				minAbsInteractionZCovariateCount = Double.parseDouble(commandLine.getOptionValue("ciz"));
			} catch (NumberFormatException ex) {
				System.out.println("Cannot not parse --covariateInteractionZ as double: " + commandLine.getOptionValue("ciz"));
				System.exit(1);
				return;
			}

			try {
				minAbsReplicationInteractionZCovariateCount = Double.parseDouble(commandLine.getOptionValue("criz"));
			} catch (NumberFormatException ex) {
				System.out.println("Cannot not parse --covariateReplicationInteractionZ as double: " + commandLine.getOptionValue("criz"));
				System.exit(1);
				return;
			}

			if (commandLine.hasOption("c")) {
				covariatesToIncludeFile = new File(commandLine.getOptionValue("c"));
			} else {
				covariatesToIncludeFile = null;
			}

			if (commandLine.hasOption("g")) {
				genesToIncludeFile = new File(commandLine.getOptionValue("g"));
			} else {
				genesToIncludeFile = null;
			}
			
			if (commandLine.hasOption("cr")) {
				covariatesReplicationToIncludeFile = new File(commandLine.getOptionValue("cr"));
			} else {
				covariatesReplicationToIncludeFile = null;
			}

			if (commandLine.hasOption("gr")) {
				genesReplicationToIncludeFile = new File(commandLine.getOptionValue("gr"));
			} else {
				genesReplicationToIncludeFile = null;
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
		BufferedWriter logWriter = new BufferedWriter(new FileWriter(outputPrefix + "_Log.txt"));

		writeAndOut("Software version: " + Main.VERSION, logWriter);
		writeAndOut("Input file: " + inputInteractionFile.getAbsolutePath(), logWriter);
		writeAndOut("Replication file: " + replicationInteractionFile.getAbsolutePath(), logWriter);
		writeAndOut("Output prefix: " + outputPrefix, logWriter);
		writeAndOut("Min interaction z-score: " + minAbsInteractionZ, logWriter);
		writeAndOut("Min replication interaction z-score: " + minAbsReplicationInteractionZ, logWriter);
		writeAndOut("Min interaction z-score covariate counter: " + minAbsInteractionZCovariateCount, logWriter);
		writeAndOut("Min replication interaction z-score covariate counter: " + minAbsReplicationInteractionZCovariateCount, logWriter);
		if (matchOnChrPos) {
			writeAndOut("Matching variants on chr-pos", logWriter);
		}
		if (covariatesToIncludeFile != null) {
			writeAndOut("Covariates to include: " + covariatesToIncludeFile.getAbsolutePath(), logWriter);
		}
		if (genesToIncludeFile != null) {
			writeAndOut("eQTL genes to include: " + genesToIncludeFile.getAbsolutePath(), logWriter);
		}
		if (covariatesReplicationToIncludeFile != null) {
			writeAndOut("Covariates replication to include: " + covariatesReplicationToIncludeFile.getAbsolutePath(), logWriter);
		}
		if (genesReplicationToIncludeFile != null) {
			writeAndOut("eQTL genes replication to include: " + genesReplicationToIncludeFile.getAbsolutePath(), logWriter);
		}
		
		writeAndOut("", logWriter);


		final HashSet<String> covariantsToInclude;
		if (covariatesToIncludeFile != null) {
			covariantsToInclude = new HashSet<String>();
			BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(covariatesToIncludeFile), "UTF-8"));
			String line;
			while ((line = reader.readLine()) != null) {
				covariantsToInclude.add(line.trim());
			}
			writeAndOut("Covariates included: " + covariantsToInclude.size(), logWriter);
			writeAndOut("", logWriter);
		} else {
			covariantsToInclude = null;
		}

		final HashSet<String> genesToInclude;
		if (genesToIncludeFile != null) {
			genesToInclude = new HashSet<String>();
			BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(genesToIncludeFile), "UTF-8"));
			String line;
			while ((line = reader.readLine()) != null) {
				genesToInclude.add(line.trim());
			}
			writeAndOut("eQTL genes included: " + genesToInclude.size(), logWriter);
			writeAndOut("", logWriter);
		} else {
			genesToInclude = null;
		}
		
		final HashSet<String> covariantsReplicationToInclude;
		if (covariatesReplicationToIncludeFile != null) {
			covariantsReplicationToInclude = new HashSet<String>();
			BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(covariatesReplicationToIncludeFile), "UTF-8"));
			String line;
			while ((line = reader.readLine()) != null) {
				covariantsReplicationToInclude.add(line.trim());
			}
			writeAndOut("Covariates replication included: " + covariantsReplicationToInclude.size(), logWriter);
			writeAndOut("", logWriter);
		} else {
			covariantsReplicationToInclude = null;
		}

		final HashSet<String> genesReplicationToInclude;
		if (genesReplicationToIncludeFile != null) {
			genesReplicationToInclude = new HashSet<String>();
			BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(genesReplicationToIncludeFile), "UTF-8"));
			String line;
			while ((line = reader.readLine()) != null) {
				genesReplicationToInclude.add(line.trim());
			}
			writeAndOut("eQTL genes replication included: " + genesReplicationToInclude.size(), logWriter);
			writeAndOut("", logWriter);
		} else {
			genesReplicationToInclude = null;
		}

		BinaryInteractionFile inputFile = BinaryInteractionFile.load(inputInteractionFile, true);
		BinaryInteractionFile replicationFile = BinaryInteractionFile.load(replicationInteractionFile, true);

		String[] row = new String[15];

		CSVWriter replicatedSameDirectionWriter = writeHeader(new File(outputPrefix + "_ReplicatedSameDirection.txt"), row);
		CSVWriter replicatedOppositeDirectionWriter = writeHeader(new File(outputPrefix + "_ReplicatedOppositeDirection.txt"), row);
		CSVWriter notReplicatedSameDirectionWriter = writeHeader(new File(outputPrefix + "_NotReplicatedSameDirection.txt"), row);
		CSVWriter notReplicatedOppositeDirectionWriter = writeHeader(new File(outputPrefix + "_NotReplicatedOppositeDirection.txt"), row);
		CSVWriter notInReplicationWriter = writeHeader(new File(outputPrefix + "_NotInReplication.txt"), row);

		int significant = 0;
		int notSignificant = 0;
		int notTestedInReplication = 0;
		int nanReplication = 0;
		int notSignificantReplicationSameDirection = 0;
		int notSignificantReplicationOppositeDirection = 0;
		int significantReplicationOppositeDirection = 0;
		int significantReplicationSameDirection = 0;

		int reporter = 0;

		LinkedHashMap<String, CovariateCount> covariateCounts = new LinkedHashMap<String, CovariateCount>(inputFile.getCovariateCount());
		for (String covariate : inputFile.getCovariates()) {
			covariateCounts.put(covariate, new CovariateCount());
		}

		for (BinaryInteractionVariant variant : inputFile.getVariants()) {

			String variantName = variant.getName();

			BinaryInteractionVariant replicationVariant;
			boolean swap;

			if (matchOnChrPos) {
				replicationVariant = replicationFile.getVariant(variant.getChr(), variant.getPos());
			} else {
				if (replicationFile.containsVariant(variantName)) {
					replicationVariant = replicationFile.getVariant(variantName);
				} else {
					replicationVariant = null;
				}
			}

			if (replicationVariant != null) {
				if (!(variant.getRefAllele() == replicationVariant.getRefAllele() && variant.getAltAllele() == replicationVariant.getAltAllele())
						&& !(variant.getRefAllele() == replicationVariant.getAltAllele() && variant.getAltAllele() == replicationVariant.getRefAllele())) {
					System.err.println("Allele mismatch!");
				}
				swap = variant.getAltAllele() != replicationVariant.getAltAllele();
			} else {
				swap = false;
			}

			//Do loop anyway to also count not replicated

			int[] genePointers = inputFile.getVariant(variantName).getGenePointers();
			genes:
			for (int genePointer : genePointers) {

				BinaryInteractionGene gene = inputFile.getGene(genePointer);
				if (genesToInclude != null && !genesToInclude.contains(gene.getName())) {
					continue genes;
				}

				covairates:
				for (Iterator<BinaryInteractionQueryResult> iterator = inputFile.readVariantGeneResults(variantName, gene.getName()); iterator.hasNext();) {

					BinaryInteractionQueryResult interaction = iterator.next();

					if (covariantsToInclude != null && !covariantsToInclude.contains(interaction.getCovariateName())) {
						continue covairates;
					}

					double metaInteractionZ = interaction.getInteractionZscores().getZscoreInteractionMeta();

					if (metaInteractionZ >= minAbsInteractionZ || metaInteractionZ <= -minAbsInteractionZ) {
						++significant;

						if (replicationVariant != null && replicationFile.containsInteraction(replicationVariant.getName(), gene.getName(), interaction.getCovariateName()) && (genesReplicationToInclude == null || genesReplicationToInclude.contains(interaction.getGeneName()) && (covariantsReplicationToInclude == null || covariantsToInclude.contains(interaction.getCovariateName())) )) {

							BinaryInteractionZscores replicationZscores = replicationFile.readInteractionResults(replicationVariant.getName(), gene.getName(), interaction.getCovariateName());
							double replicationInteractionZscore = replicationZscores.getZscoreInteractionMeta();

							BinaryInteractionQtlZscores replicationQtlRes = replicationFile.readQtlResults(replicationVariant.getName(), gene.getName());

							if (!Double.isNaN(replicationInteractionZscore)) {

								if (swap) {
									replicationInteractionZscore *= -1;
								}

								if (replicationInteractionZscore <= -minAbsReplicationInteractionZ || replicationInteractionZscore >= minAbsReplicationInteractionZ) {
									if (metaInteractionZ * replicationInteractionZscore >= 0) {
										++significantReplicationSameDirection;


										writeInteraction(row, variantName, gene, interaction, variant, replicationQtlRes, replicationZscores, swap, replicatedSameDirectionWriter);

									} else {
										++significantReplicationOppositeDirection;

										writeInteraction(row, variantName, gene, interaction, variant, replicationQtlRes, replicationZscores, swap, replicatedOppositeDirectionWriter);
									}
								} else {
									if (metaInteractionZ * replicationInteractionZscore >= 0) {
										++notSignificantReplicationSameDirection;
										writeInteraction(row, variantName, gene, interaction, variant, replicationQtlRes, replicationZscores, swap, notReplicatedSameDirectionWriter);
									} else {
										++notSignificantReplicationOppositeDirection;
										writeInteraction(row, variantName, gene, interaction, variant, replicationQtlRes, replicationZscores, swap, notReplicatedOppositeDirectionWriter);
									}
								}
							} else {
								writeInteraction(row, variantName, gene, interaction, variant, replicationQtlRes, replicationZscores, swap, notInReplicationWriter);
								++nanReplication;
							}

						} else {
							writeInteraction(row, variantName, gene, interaction, variant, null, null, swap, notInReplicationWriter);
							++notTestedInReplication;
						}

					} else {
						++notSignificant;
					}

					if (metaInteractionZ >= minAbsInteractionZCovariateCount || metaInteractionZ <= -minAbsInteractionZCovariateCount) {

						CovariateCount thisCovariateCounts = covariateCounts.get(interaction.getCovariateName());
						thisCovariateCounts.incrementCovariateSignificant();

						if (replicationVariant != null && replicationFile.containsInteraction(replicationVariant.getName(), gene.getName(), interaction.getCovariateName())) {

							BinaryInteractionZscores replicationZscores = replicationFile.readInteractionResults(replicationVariant.getName(), gene.getName(), interaction.getCovariateName());
							double replicationInteractionZscore = replicationZscores.getZscoreInteractionMeta();

							if (!Double.isNaN(replicationInteractionZscore)) {

								if (swap) {
									replicationInteractionZscore *= -1;
								}

								if (replicationInteractionZscore <= -minAbsReplicationInteractionZCovariateCount || replicationInteractionZscore >= minAbsReplicationInteractionZCovariateCount) {
									if (metaInteractionZ * replicationInteractionZscore >= 0) {
										thisCovariateCounts.incrementReplicatedSameDirection();

									} else {
										thisCovariateCounts.incrementReplicatedOppositeDirection();
									}
								} else {
									if (metaInteractionZ * replicationInteractionZscore >= 0) {
										thisCovariateCounts.incrementNotReplicatedSameDirection();
									} else {
										thisCovariateCounts.incrementNotReplicatedOppositeDirection();
									}
								}

							} else {
							}

						} else {
						}

					} else {
					}


				}

				++reporter;
				if (reporter % 500 == 0) {
					System.out.println("Parsed " + reporter + " of " + inputFile.getVariantGeneCombinations() + " variant-gene combinations");
				}

			}

		}

		replicatedSameDirectionWriter.close();
		replicatedOppositeDirectionWriter.close();
		notReplicatedSameDirectionWriter.close();
		notReplicatedOppositeDirectionWriter.close();
		notInReplicationWriter.close();

		writeCovaraiteCounts(new File(outputPrefix + "_CovariateCounts.txt"), covariateCounts);

		NumberFormat numberFormat = NumberFormat.getInstance();
		numberFormat.setMinimumFractionDigits(0);
		numberFormat.setMaximumFractionDigits(2);

		writeAndOut("", logWriter);
		writeAndOut("Total number of interactions: " + numberFormat.format(notSignificant + significant), logWriter);
		writeAndOut(" - Not significant: " + numberFormat.format(notSignificant) + " (" + numberFormat.format(notSignificant * 100d / (notSignificant + significant)) + "%)", logWriter);
		writeAndOut(" - Significant: " + numberFormat.format(significant) + " (" + numberFormat.format(significant * 100d / (notSignificant + significant)) + "%)", logWriter);
		writeAndOut("  * Not in replication: " + numberFormat.format(notTestedInReplication) + " (" + numberFormat.format(notTestedInReplication * 100d / significant) + "%)", logWriter);
		writeAndOut("  * NaN in replication: " + numberFormat.format(nanReplication) + " (" + numberFormat.format(nanReplication * 100d / significant) + "%)", logWriter);
		writeAndOut("  * Not significant in replication: " + numberFormat.format(notSignificantReplicationSameDirection + notSignificantReplicationOppositeDirection) + " (" + numberFormat.format((notSignificantReplicationSameDirection + notSignificantReplicationOppositeDirection) * 100d / significant) + "%)", logWriter);
		writeAndOut("   # Same direction: " + numberFormat.format(notSignificantReplicationSameDirection) + " (" + numberFormat.format(notSignificantReplicationSameDirection * 100d / (notSignificantReplicationSameDirection + notSignificantReplicationOppositeDirection)) + "%)", logWriter);
		writeAndOut("   # Opposite direction: " + numberFormat.format(notSignificantReplicationOppositeDirection) + " (" + numberFormat.format(notSignificantReplicationOppositeDirection * 100d / (notSignificantReplicationSameDirection + notSignificantReplicationOppositeDirection)) + "%)", logWriter);
		writeAndOut("  * Significant in replication: " + numberFormat.format(significantReplicationSameDirection + significantReplicationOppositeDirection) + " (" + numberFormat.format((significantReplicationSameDirection + significantReplicationOppositeDirection) * 100d / significant) + "%)", logWriter);
		writeAndOut("   # Same direction: " + numberFormat.format(significantReplicationSameDirection) + " (" + numberFormat.format(significantReplicationSameDirection * 100d / (significantReplicationSameDirection + significantReplicationOppositeDirection)) + "%)", logWriter);
		writeAndOut("   # Opposite direction: " + numberFormat.format(significantReplicationOppositeDirection) + " (" + numberFormat.format(significantReplicationOppositeDirection * 100d / (significantReplicationSameDirection + significantReplicationOppositeDirection)) + "%)", logWriter);

		logWriter.close();

	}

	private static void writeInteraction(String[] row, String variantName, BinaryInteractionGene gene, BinaryInteractionQueryResult interation, BinaryInteractionVariant variant, BinaryInteractionQtlZscores replicationQtlRes, BinaryInteractionZscores replicationZscores, boolean swap, CSVWriter interactionWriter) {
		int c = 0;
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
		row[c++] = replicationQtlRes == null ? "NaN" : String.valueOf(replicationQtlRes.getMetaZscore() * (swap ? -1 : 1));
		row[c++] = replicationZscores == null ? "NaN" : String.valueOf(replicationZscores.getZscoreSnpMeta() * (swap ? -1 : 1));
		row[c++] = replicationZscores == null ? "NaN" : String.valueOf(replicationZscores.getZscoreCovariateMeta());
		row[c++] = replicationZscores == null ? "NaN" : String.valueOf(replicationZscores.getZscoreInteractionMeta() * (swap ? -1 : 1));
		interactionWriter.writeNext(row);
	}

	private static CSVWriter writeHeader(File file, String[] row) throws IOException {
		CSVWriter replicatedSameDirectionWriter = new CSVWriter(new BufferedWriter(new FileWriter(file)), '\t', '\0', '\0');
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
		return replicatedSameDirectionWriter;
	}

	private static void writeCovaraiteCounts(File file, LinkedHashMap<String, CovariateCount> covariateCounts) throws IOException {

		CSVWriter covariateCountWriter = new CSVWriter(new BufferedWriter(new FileWriter(file)), '\t', '\0', '\0');
		int c = 0;
		String[] row2 = new String[6];
		row2[c++] = "Covariate";
		row2[c++] = "Significant";
		row2[c++] = "ReplicatedSameDirection";
		row2[c++] = "ReplicatedOppositeDirection";
		row2[c++] = "NotReplicateSameDirection";
		row2[c++] = "NotReplicatedOppositeDirection";
		covariateCountWriter.writeNext(row2);

		for (Map.Entry<String, CovariateCount> covariateEntry : covariateCounts.entrySet()) {

			CovariateCount thisCounts = covariateEntry.getValue();

			c = 0;
			row2[c++] = covariateEntry.getKey();
			row2[c++] = String.valueOf(thisCounts.getCovariateSignificant());
			row2[c++] = String.valueOf(thisCounts.getReplicatedSameDirection());
			row2[c++] = String.valueOf(thisCounts.getReplicatedOppositeDirection());
			row2[c++] = String.valueOf(thisCounts.getNotReplicatedSameDirection());
			row2[c++] = String.valueOf(thisCounts.getNotReplicatedOppositeDirection());
			covariateCountWriter.writeNext(row2);

		}

		covariateCountWriter.close();

	}

	private static void writeAndOut(String message, Writer writer) throws IOException {
		writer.append(message);
		writer.append('\n');
		System.out.println(message);
	}

	private static class CovariateCount {

		private int covariateSignificant = 0;
		private int replicatedSameDirection = 0;
		private int replicatedOppositeDirection = 0;
		private int notReplicatedSameDirection = 0;
		private int notReplicatedOppositeDirection = 0;

		public int getCovariateSignificant() {
			return covariateSignificant;
		}

		public int getReplicatedSameDirection() {
			return replicatedSameDirection;
		}

		public int getReplicatedOppositeDirection() {
			return replicatedOppositeDirection;
		}

		public int getNotReplicatedSameDirection() {
			return notReplicatedSameDirection;
		}

		public int getNotReplicatedOppositeDirection() {
			return notReplicatedOppositeDirection;
		}

		public void incrementCovariateSignificant() {
			covariateSignificant++;
		}

		public void incrementReplicatedSameDirection() {
			replicatedSameDirection++;
		}

		public void incrementReplicatedOppositeDirection() {
			replicatedOppositeDirection++;
		}

		public void incrementNotReplicatedSameDirection() {
			notReplicatedSameDirection++;
		}

		public void incrementNotReplicatedOppositeDirection() {
			notReplicatedOppositeDirection++;
		}
	}
}
