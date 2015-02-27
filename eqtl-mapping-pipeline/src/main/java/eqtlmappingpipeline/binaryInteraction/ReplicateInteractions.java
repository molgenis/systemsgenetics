package eqtlmappingpipeline.binaryInteraction;

import java.io.File;
import java.io.FileNotFoundException;
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

		try {
			final CommandLine commandLine = new PosixParser().parse(OPTIONS, args, false);

			inputInteractionFile = new File(commandLine.getOptionValue("i"));
			replicationInteractionFile = new File(commandLine.getOptionValue("r"));

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
		System.out.println("Min interaction z-score: " + minAbsInteractionZ);
		System.out.println("Min replication interaction z-score: " + minAbsReplicationInteractionZ);
		if(matchOnChrPos){
			System.out.println("Matching variants on chr-pos");
		}
		System.out.println("");

		BinaryInteractionFile inputFile = BinaryInteractionFile.load(inputInteractionFile, true);
		BinaryInteractionFile replicationFile = BinaryInteractionFile.load(replicationInteractionFile, true);

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
			
			if(matchOnChrPos){
				replicationVariant = replicationFile.getVariant(variant.getChr(), variant.getPos());
			} else {
				if(replicationFile.containsVariant(variantName)){
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
	
							if(variant.getAltAllele() != replicationVariant.getAltAllele()){
								replicationInteractionZscore *= -1;
							}
							
							if(replicationInteractionZscore <= -minAbsReplicationInteractionZ || replicationInteractionZscore >= minAbsReplicationInteractionZ){
								if(metaInteractionZ * replicationInteractionZscore >= 0){
									++significantReplicationSameDirection;
								} else {
									++significantReplicationOppositeDirection;
								}
							} else {
								if(metaInteractionZ * replicationInteractionZscore >= 0){
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
				if(reporter % 500 == 0){
					System.out.println("Parsed " + reporter + " of " + inputFile.getVariantGeneCombinations() + " variant-gene combinations");
				}
					
			}
		}
		
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
