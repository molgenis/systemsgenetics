package eqtlmappingpipeline.binaryInteraction;

import au.com.bytecode.opencsv.CSVWriter;
import gnu.trove.map.hash.TObjectDoubleHashMap;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashSet;
import java.util.Iterator;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import umcg.genetica.io.binInteraction.BinaryInteractionFile;
import umcg.genetica.io.binInteraction.BinaryInteractionFileException;
import umcg.genetica.io.binInteraction.BinaryInteractionQueryResult;
import umcg.genetica.io.binInteraction.gene.BinaryInteractionGene;
import umcg.genetica.io.binInteraction.variant.BinaryInteractionVariant;

/**
 *
 * @author Patrick Deelen
 */
public class CovariateImportance {

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
		OptionBuilder.withDescription("Ouput file");
		OptionBuilder.withLongOpt("output");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("o"));

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

	}

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws FileNotFoundException, IOException, BinaryInteractionFileException {

		final File inputInteractionFile;
		final File outputFile;
		final File covariatesToIncludeFile;
		final File genesToIncludeFile;

		try {
			final CommandLine commandLine = new PosixParser().parse(OPTIONS, args, false);

			inputInteractionFile = new File(commandLine.getOptionValue("i"));
			outputFile = new File(commandLine.getOptionValue("o"));

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

		} catch (ParseException ex) {
			System.err.println("Invalid command line arguments: ");
			System.err.println(ex.getMessage());
			System.err.println();
			new HelpFormatter().printHelp(" ", OPTIONS);
			System.exit(1);
			return;
		}

		System.out.println("Input file: " + inputInteractionFile.getAbsolutePath());
		System.out.println("Output file: " + outputFile);
		if (covariatesToIncludeFile != null) {
			System.out.println("Covariates to include: " + covariatesToIncludeFile.getAbsolutePath());
		}
		if (genesToIncludeFile != null) {
			System.out.println("eQTL genes to include: " + genesToIncludeFile.getAbsolutePath());
		}
		
		final HashSet<String> genesToInclude;
		if (genesToIncludeFile != null) {
			genesToInclude = new HashSet<String>();
			BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(genesToIncludeFile), "UTF-8"));
			String line;
			while ((line = reader.readLine()) != null) {
				genesToInclude.add(line.trim());
			}
			System.out.println("eQTL genes included: " + genesToInclude.size());
			System.out.println("");
		} else {
			genesToInclude = null;
		}

		final HashSet<String> covariantsToInclude;
		if (covariatesToIncludeFile != null) {
			covariantsToInclude = new HashSet<String>();
			BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(covariatesToIncludeFile), "UTF-8"));
			String line;
			while ((line = reader.readLine()) != null) {
				covariantsToInclude.add(line.trim());
			}
			System.out.println("Covariates included: " + covariantsToInclude.size());
			System.out.println();
		} else {
			covariantsToInclude = null;
		}

		BinaryInteractionFile inputFile = BinaryInteractionFile.load(inputInteractionFile, true);

		TObjectDoubleHashMap<String> sumChi2 = new TObjectDoubleHashMap<String>(20000, 0.75f, Double.NaN);

		int reporter = 0;
		for (BinaryInteractionVariant variant : inputFile.getVariants()) {

			String variantName = variant.getName();
			int[] genePointers = variant.getGenePointers();

			genes:
			for (int genePointer : genePointers) {

				BinaryInteractionGene gene = inputFile.getGene(genePointer);
				
				if (genesToInclude != null && !genesToInclude.contains(gene.getName())) {
					continue genes;
				}

				covariates:
				for (Iterator<BinaryInteractionQueryResult> iterator = inputFile.readVariantGeneResults(variantName, gene.getName()); iterator.hasNext();) {

					BinaryInteractionQueryResult interation = iterator.next();

					if (covariantsToInclude != null && !covariantsToInclude.contains(interation.getCovariateName())) {
						continue covariates;
					}

					double metaZ = interation.getInteractionZscores().getZscoreInteractionMeta();
					if (Double.isNaN(metaZ)) {
						continue covariates;
					}
					double chi2 = metaZ * metaZ;
					sumChi2.adjustOrPutValue(interation.getCovariateName(), chi2, chi2);

				}

				++reporter;
				if (reporter % 500 == 0) {
					System.out.println("Parsed " + reporter + " of " + inputFile.getVariantGeneCombinations() + " variant-gene combinations");
				}

			}



		}

		CSVWriter outputWriter = new CSVWriter(new BufferedWriter(new FileWriter(outputFile)), '\t', '\0', '\0');
		String[] row = new String[2];
		int c = 0;
		row[c++] = "Covariate";
		row[c++] = "sumChi2";
		outputWriter.writeNext(row);

		for (String covariate : inputFile.getCovariates()) {
			c = 0;
			row[c++] = covariate;
			row[c++] = String.valueOf(sumChi2.get(covariate));
			outputWriter.writeNext(row);
		}

		outputWriter.close();
		System.out.println("Done");

	}
}