package eqtlmappingpipeline.binaryInteraction;

import au.com.bytecode.opencsv.CSVWriter;
import eqtlmappingpipeline.Main;
import gnu.trove.map.hash.TObjectIntHashMap;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Writer;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
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
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author Patrick Deelen
 */
public class InvestigateCovariate {

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

		OptionBuilder.withArgName("string");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Covariate name");
		OptionBuilder.withLongOpt("queryCovariate");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("qc"));

		OptionBuilder.withDescription("If set match variant on chr-pos");
		OptionBuilder.withLongOpt("chrPos");
		OPTIONS.addOption(OptionBuilder.create("cp"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("File with covariates to include in analysis");
		OptionBuilder.withLongOpt("covariats");
		OPTIONS.addOption(OptionBuilder.create("c"));

	}

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException, FileNotFoundException, BinaryInteractionFileException {

		final File inputInteractionFile;
		final File replicationInteractionFile;
		final double minAbsInteractionZ;
		final double minAbsReplicationInteractionZ;
		final boolean matchOnChrPos;
		final String outputPrefix;
		final String queryCovariateName;
		final File covariatesToIncludeFile;

		try {
			final CommandLine commandLine = new PosixParser().parse(OPTIONS, args, false);

			inputInteractionFile = new File(commandLine.getOptionValue("i"));
			replicationInteractionFile = new File(commandLine.getOptionValue("r"));
			outputPrefix = commandLine.getOptionValue("o");
			queryCovariateName = commandLine.getOptionValue("qc");

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

			if (commandLine.hasOption("c")) {
				covariatesToIncludeFile = new File(commandLine.getOptionValue("c"));
			} else {
				covariatesToIncludeFile = null;
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
		writeAndOut("Query covariate: " + queryCovariateName, logWriter);
		writeAndOut("Output prefix: " + outputPrefix, logWriter);
		writeAndOut("Min interaction z-score: " + minAbsInteractionZ, logWriter);
		writeAndOut("Min replication interaction z-score: " + minAbsReplicationInteractionZ, logWriter);
		if (matchOnChrPos) {
			writeAndOut("Matching variants on chr-pos", logWriter);
		}
		if (covariatesToIncludeFile != null) {
			writeAndOut("Covariates to include: " + covariatesToIncludeFile.getAbsolutePath(), logWriter);
		}
		writeAndOut("", logWriter);

		final HashSet<String> covariantsToIncluded;
		if (covariatesToIncludeFile != null) {
			covariantsToIncluded = new HashSet<String>();
			BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(covariatesToIncludeFile), "UTF-8"));
			String line;
			while ((line = reader.readLine()) != null) {
				covariantsToIncluded.add(line.trim());
			}
			writeAndOut("Covariates included: " + covariantsToIncluded.size(), logWriter);
			writeAndOut("", logWriter);

			if (!covariantsToIncluded.contains(queryCovariateName)) {
				System.err.println("Query covariate not in include list");
				System.exit(1);
				return;
			}

		} else {
			covariantsToIncluded = null;
		}

		BinaryInteractionFile inputFile = BinaryInteractionFile.load(inputInteractionFile, true);
		BinaryInteractionFile replicationFile = BinaryInteractionFile.load(replicationInteractionFile, true);

		LinkedHashSet<String> genesOfInterest = new LinkedHashSet<String>();

		if (!inputFile.containsCovariant(queryCovariateName)) {
			System.err.println("Covariate not found in input data");
			System.exit(1);
			return;
		}

		if (!replicationFile.containsCovariant(queryCovariateName)) {
			System.err.println("Covariate not found in replication data");
			System.exit(1);
			return;
		}

		variants:
		for (final BinaryInteractionVariant variant : inputFile.getVariants()) {

			final String variantName = variant.getName();

			final BinaryInteractionVariant replicationVariant;

			if (matchOnChrPos) {
				replicationVariant = replicationFile.getVariant(variant.getChr(), variant.getPos());
				if (replicationVariant == null) {
					continue variants;
				}
			} else {
				if (replicationFile.containsVariant(variantName)) {
					replicationVariant = replicationFile.getVariant(variantName);
				} else {
					continue variants;
				}
			}

			//Only do if replication variant has been found

			if (!(variant.getRefAllele() == replicationVariant.getRefAllele() && variant.getAltAllele() == replicationVariant.getAltAllele())
					&& !(variant.getRefAllele() == replicationVariant.getAltAllele() && variant.getAltAllele() == replicationVariant.getRefAllele())) {
				System.err.println("Allele mismatch!");
			}
			final boolean swap = variant.getAltAllele() != replicationVariant.getAltAllele();

			final int[] genePointers = inputFile.getVariant(variantName).getGenePointers();

			genes:
			for (int genePointer : genePointers) {

				final BinaryInteractionGene gene = inputFile.getGene(genePointer);

				if (!inputFile.containsInteraction(variantName, gene.getName(), queryCovariateName)) {
					continue genes;
				}

				if (!replicationFile.containsInteraction(replicationVariant.getName(), gene.getName(), queryCovariateName)) {
					continue genes;
				}

				if (!replicationFile.containsInteraction(replicationVariant.getName(), gene.getName(), queryCovariateName)) {
					continue genes;
				}

				final BinaryInteractionZscores inputInteractionResult = inputFile.readInteractionResults(variantName, gene.getName(), queryCovariateName);
				final double inputInteractionZ = inputInteractionResult.getZscoreInteractionMeta();

				if (Double.isNaN(inputInteractionZ)) {
					continue genes;
				}

				if (!(inputInteractionZ <= -minAbsInteractionZ || inputInteractionZ >= minAbsInteractionZ)) {
					continue genes;
				}
				
				if(!replicationFile.containsInteraction(replicationVariant.getName(), gene.getName(), queryCovariateName)){
					continue genes;
				}

				final BinaryInteractionZscores replicationInteractionResult = replicationFile.readInteractionResults(replicationVariant.getName(), gene.getName(), queryCovariateName);
				double replicationInteractionZ = replicationInteractionResult.getZscoreInteractionMeta();

				if (Double.isNaN(replicationInteractionZ)) {
					continue genes;
				}

				if (!(replicationInteractionZ <= -minAbsReplicationInteractionZ || replicationInteractionZ >= minAbsReplicationInteractionZ)) {
					continue genes;
				}

				//If here then discovery and replication significant

				if (swap) {
					replicationInteractionZ *= -1;
				}

				if (inputInteractionZ * replicationInteractionZ >= 0) {
					//Same direction
					genesOfInterest.add(gene.getName());
				}



			}

		}
		
		System.out.println("Number of genes of interest: " + genesOfInterest.size());

		TObjectIntHashMap<String> covaraitesOfInterestCount = new TObjectIntHashMap<String>();
		TObjectIntHashMap<String> genesOfInterestCount = new TObjectIntHashMap<String>();
		LinkedHashSet<String> covaraitesOfInterest = new LinkedHashSet<String>();

		//Here we now know which genes are of interest.
		//We are now going to search for other covariates that are significant for any of these genes
		for (String geneName : genesOfInterest) {
			BinaryInteractionGene gene = inputFile.getGene(geneName);

			variants:
			for (int variantPointer : gene.getVariantPointers()) {
				BinaryInteractionVariant variant = inputFile.getVariant(variantPointer);

				final String variantName = variant.getName();

				final BinaryInteractionVariant replicationVariant;

				if (matchOnChrPos) {
					replicationVariant = replicationFile.getVariant(variant.getChr(), variant.getPos());
					if (replicationVariant == null) {
						continue variants;
					}
				} else {
					if (replicationFile.containsVariant(variantName)) {
						replicationVariant = replicationFile.getVariant(variantName);
					} else {
						continue variants;
					}
				}

				//Only do if replication variant has been found

				if (!(variant.getRefAllele() == replicationVariant.getRefAllele() && variant.getAltAllele() == replicationVariant.getAltAllele())
						&& !(variant.getRefAllele() == replicationVariant.getAltAllele() && variant.getAltAllele() == replicationVariant.getRefAllele())) {
					System.err.println("Allele mismatch!");
				}
				final boolean swap = variant.getAltAllele() != replicationVariant.getAltAllele();

				covairates:
				for (Iterator<BinaryInteractionQueryResult> iterator = inputFile.readVariantGeneResults(variantName, gene.getName()); iterator.hasNext();) {

					BinaryInteractionQueryResult interaction = iterator.next();

					if (covariantsToIncluded != null && !covariantsToIncluded.contains(interaction.getCovariateName())) {
						continue covairates;
					}

					final BinaryInteractionZscores inputInteractionResult = interaction.getInteractionZscores();
					final double inputInteractionZ = inputInteractionResult.getZscoreInteractionMeta();

					if (Double.isNaN(inputInteractionZ)) {
						continue;
					}

					if (!(inputInteractionZ <= -minAbsInteractionZ || inputInteractionZ >= minAbsInteractionZ)) {
						continue covairates;
					}
					
					if(!replicationFile.containsInteraction(replicationVariant.getName(), gene.getName(), interaction.getCovariateName())){
						continue covairates;
					}

					final BinaryInteractionZscores replicationInteractionResult = replicationFile.readInteractionResults(replicationVariant.getName(), gene.getName(), interaction.getCovariateName());
					double replicationInteractionZ = replicationInteractionResult.getZscoreInteractionMeta();

					if (Double.isNaN(replicationInteractionZ)) {
						continue covairates;
					}

					if (!(replicationInteractionZ <= -minAbsReplicationInteractionZ || replicationInteractionZ >= minAbsReplicationInteractionZ)) {
						continue covairates;
					}

					//If here then discovery and replication significant

					if (swap) {
						replicationInteractionZ *= -1;
					}

					if (inputInteractionZ * replicationInteractionZ >= 0) {
						//Same direction
						covaraitesOfInterestCount.adjustOrPutValue(interaction.getCovariateName(), 1, 1);
						covaraitesOfInterest.add(interaction.getCovariateName());
						genesOfInterestCount.adjustOrPutValue(geneName, 1, 1);
					}

				}

			}

		}
		//We now also know which other covariates are of interest
		
		System.out.println("Number of covariates of interest: " + covaraitesOfInterest.size());
		
		
		writeCounts(genesOfInterest, genesOfInterestCount, new File(outputPrefix + "_Genes.txt"));
		writeCounts(covaraitesOfInterest, covaraitesOfInterestCount, new File(outputPrefix + "_Covariates.txt"));

		DoubleMatrixDataset<String, String> interactionZscores = new DoubleMatrixDataset<String, String>(covaraitesOfInterest, genesOfInterest);

		DoubleMatrixDataset<String, String> replicationInteractionZscores = new DoubleMatrixDataset<String, String>(covaraitesOfInterest, genesOfInterest);

		for (String geneName : genesOfInterest) {

			BinaryInteractionGene gene = inputFile.getGene(geneName);

			variants:
			for (int variantPointer : gene.getVariantPointers()) {
				BinaryInteractionVariant variant = inputFile.getVariant(variantPointer);

				final String variantName = variant.getName();

				final BinaryInteractionVariant replicationVariant;

				if (matchOnChrPos) {
					replicationVariant = replicationFile.getVariant(variant.getChr(), variant.getPos());
					if (replicationVariant == null) {
						continue variants;
					}
				} else {
					if (replicationFile.containsVariant(variantName)) {
						replicationVariant = replicationFile.getVariant(variantName);
					} else {
						continue variants;
					}
				}

				//Only do if replication variant has been found

				if (!(variant.getRefAllele() == replicationVariant.getRefAllele() && variant.getAltAllele() == replicationVariant.getAltAllele())
						&& !(variant.getRefAllele() == replicationVariant.getAltAllele() && variant.getAltAllele() == replicationVariant.getRefAllele())) {
					System.err.println("Allele mismatch!");
				}
				final boolean swap = variant.getAltAllele() != replicationVariant.getAltAllele();


				for (String covariateName : covaraitesOfInterest) {

					final BinaryInteractionZscores inputInteractionResult = inputFile.readInteractionResults(variantName, geneName, covariateName);

					//System.out.println(covariateName +  "-"  + geneName +  "-"  + inputInteractionResult.getZscoreInteractionMeta());
					interactionZscores.setElement(covariateName, geneName, inputInteractionResult.getZscoreInteractionMeta());

					final BinaryInteractionZscores replicationInteractionResult = replicationFile.readInteractionResults(replicationVariant.getName(), gene.getName(), covariateName);
					double replicationInteractionZ = replicationInteractionResult.getZscoreInteractionMeta();
					if (swap) {
						replicationInteractionZ *= -1;
					}
					replicationInteractionZscores.setElement(covariateName, geneName, replicationInteractionZ);

				}


			}
		}

		interactionZscores.save(outputPrefix + "_InteractionMatrix.txt");
		replicationInteractionZscores.save(outputPrefix + "_ReplicationInteractionMatrix.txt");
		
		logWriter.close();

	}

	private static void writeCounts(LinkedHashSet<String> elements, TObjectIntHashMap<String> counts, File file) throws IOException {
		CSVWriter writer = new CSVWriter(new BufferedWriter(new FileWriter(file)), '\t', '\0', '\0');

		String[] row = new String[2];

		for (String elementName : elements) {
			row[0] = elementName;
			row[1] = String.valueOf(counts.get(elementName));
			writer.writeNext(row);
		}
		
		writer.close();

	}

	private static void writeAndOut(String message, Writer writer) throws IOException {
		writer.append(message);
		writer.append('\n');
		System.out.println(message);
	}
}
