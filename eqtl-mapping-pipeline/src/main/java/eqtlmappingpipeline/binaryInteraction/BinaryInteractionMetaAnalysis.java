package eqtlmappingpipeline.binaryInteraction;

import eqtlmappingpipeline.Main;
import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import umcg.genetica.io.binInteraction.BinaryInteractionCohort;
import umcg.genetica.io.binInteraction.BinaryInteractionFile;
import umcg.genetica.io.binInteraction.BinaryInteractionFileCreator;
import umcg.genetica.io.binInteraction.BinaryInteractionQtlZscores;
import umcg.genetica.io.binInteraction.BinaryInteractionZscores;
import umcg.genetica.io.binInteraction.gene.BinaryInteractionGene;
import umcg.genetica.io.binInteraction.gene.BinaryInteractionGeneCreator;
import umcg.genetica.io.binInteraction.variant.BinaryInteractionVariant;
import umcg.genetica.io.binInteraction.variant.BinaryInteractionVariantCreator;

/**
 *
 * @author Patrick Deelen
 */
public class BinaryInteractionMetaAnalysis {

	private static final String VERSION = Main.VERSION;
	private static final String HEADER =
			"  /---------------------------------------\\\n"
			+ "  |   Binary interaction meta analysis    |\n"
			+ "  |                                       |\n"
			+ "  |             Patrick Deelen            |\n"
			+ "  |        patrickdeelen@gmail.com        |\n"
			+ "  |                                       |\n"
			+ "  |   Dasha Zhernakova, Marc Jan Bonder   |\n"
			+ "  |      Lude Franke, Morris Swertz       |\n"
			+ "  |                                       |\n"
			+ "  |     Genomics Coordication Center      |\n"
			+ "  |        Department of Genetics         |\n"
			+ "  |  University Medical Center Groningen  |\n"
			+ "  \\---------------------------------------/";
	private static final DateFormat DATE_TIME_FORMAT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
	private static final Date currentDataTime = new Date();
	private static final Options OPTIONS;

	static {

		OPTIONS = new Options();

		Option option;

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Binary interaction file");
		OptionBuilder.withLongOpt("input");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("i"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Output file");
		OptionBuilder.withLongOpt("output");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("o"));

	}

	public static void main(String[] args) throws UnsupportedEncodingException, IOException, Exception {

		System.out.println(HEADER);
		System.out.println();
		System.out.println("     --- Version: " + VERSION + " ---");
		System.out.println();

		System.out.println("Current date and time: " + DATE_TIME_FORMAT.format(currentDataTime));
		System.out.println();

		System.out.flush(); //flush to make sure header is before errors
		try {
			Thread.sleep(25); //Allows flush to complete
		} catch (InterruptedException ex) {
		}

		final File[] inputInteractionFiles;
		final File outputFile;

		try {
			final CommandLine commandLine = new PosixParser().parse(OPTIONS, args, false);

			String[] inputPaths = commandLine.getOptionValues("i");
			inputInteractionFiles = new File[inputPaths.length];

			for (int i = 0; i < inputPaths.length; ++i) {
				inputInteractionFiles[i] = new File(inputPaths[i]);
			}

			outputFile = new File(commandLine.getOptionValue("o"));


		} catch (ParseException ex) {
			System.err.println("Invalid command line arguments: ");
			System.err.println(ex.getMessage());
			System.err.println();
			new HelpFormatter().printHelp(" ", OPTIONS);
			System.exit(1);
			return;
		}

		if (inputInteractionFiles.length < 1) {
			System.err.println("Supply atleast two input files for meta analysis");
			System.exit(1);
			return;
		}

		System.out.println("Input files: ");
		for (File file : inputInteractionFiles) {
			System.out.println(" * " + file.getAbsolutePath());
		}
		System.out.println("Output file: " + outputFile);

		BinaryInteractionFile[] binaryInteractionFiles = new BinaryInteractionFile[inputInteractionFiles.length];
		for (int i = 0; i < inputInteractionFiles.length; ++i) {
			binaryInteractionFiles[i] = BinaryInteractionFile.load(inputInteractionFiles[i]);
		}

		LinkedHashMap<String, BinaryInteractionVariantCreator> variants = new LinkedHashMap<String, BinaryInteractionVariantCreator>();
		LinkedHashMap<String, BinaryInteractionGeneCreator> genes = new LinkedHashMap<String, BinaryInteractionGeneCreator>();
		LinkedHashSet<String> covariates = new LinkedHashSet<String>();
		LinkedHashMap<String, BinaryInteractionCohort> cohorts = new LinkedHashMap<String, BinaryInteractionCohort>();
		LinkedHashSet<VariantGene> variantGenes = new LinkedHashSet<VariantGene>();

		boolean allStoreQtl = true;
		boolean allStoreFlipped = true;
		boolean allInteractions = true;

		for (BinaryInteractionFile binaryInteractionFile : binaryInteractionFiles) {

			if (!binaryInteractionFile.isNormalQtlStored()) {
				allStoreQtl = false;
			}

			if (!binaryInteractionFile.isFlippedZscoreStored()) {
				allStoreFlipped = false;
			}

			if (!binaryInteractionFile.areAllCovariatesTestedForAllVariantGenes()) {
				allInteractions = false;
			}

			List<BinaryInteractionGene> fileGenes = binaryInteractionFile.getGenes();

			for (BinaryInteractionVariant variant : binaryInteractionFile.getVariants()) {
				if (!variants.containsKey(variant.getName())) {
					variants.put(variant.getName(), new BinaryInteractionVariantCreator(variant.getName(), variant.getChr(), variant.getPos(), variant.getRefAllele(), variant.getAltAllele()));
				}
				for (int geneIndex : variant.getGenePointers()) {
					variantGenes.add(new VariantGene(variant.getName(), fileGenes.get(geneIndex).getName()));
				}
			}

			for (BinaryInteractionGene gene : fileGenes) {
				if (!genes.containsKey(gene.getName())) {
					genes.put(gene.getName(), new BinaryInteractionGeneCreator(gene.getName(), gene.getChr(), gene.getStart(), gene.getEnd()));
				}
			}

			for (String covariate : binaryInteractionFile.getCovariates()) {
				covariates.add(covariate);
			}

			for (BinaryInteractionCohort cohort : binaryInteractionFile.getCohorts()) {
				String cohortName = cohort.getName();
				if (cohorts.containsKey(cohortName)) {
					int i = 0;
					while (cohorts.containsKey(cohortName)) {
						cohortName = cohort.getName() + "_" + i++;
					}
					System.out.println("WARNING: cohort " + cohort.getName() + " found multiple times. Renaming to " + cohortName + " for file: " + binaryInteractionFile.getInteractionFile().getAbsolutePath());
				}
			}

		}

		final int cohortCount = cohorts.size();

		if (!allInteractions) {
			System.out.println("One or more of the input files do not contain results for all covariates for all variant-gene combinations");
		}

		if (!allStoreFlipped) {
			System.out.println("One or more of the input files do not store filled z-score. These will now also not be stored in the output file");
		}

		if (!allStoreQtl) {
			System.out.println("One or more of the input files do not normal qtl effects. These will now also not be stored in the output file");
		}

		BinaryInteractionFileCreator outputCreator = new BinaryInteractionFileCreator(
				outputFile,
				variants.values().toArray(new BinaryInteractionVariantCreator[variants.size()]),
				genes.values().toArray(new BinaryInteractionGeneCreator[genes.size()]),
				cohorts.values().toArray(new BinaryInteractionCohort[cohorts.size()]),
				covariates.toArray(new String[covariates.size()]),
				true,
				true,
				allStoreQtl,
				allStoreFlipped);

		for (VariantGene variantGene : variantGenes) {
			outputCreator.addTestedVariantGene(variantGene.getVariantName(), variantGene.getGeneName());
		}

		BinaryInteractionFile output = outputCreator.create();

		for (String variant : variants.keySet()) {

			for (String gene : genes.keySet()) {

				if (allStoreQtl) {

					int[] sampleCountsQtl = new int[cohortCount];
					double[] zscoresQtl = new double[cohortCount];

					int i = 0;
					for (BinaryInteractionFile binaryInteractionFile : binaryInteractionFiles) {

						if (binaryInteractionFile.containsVariantGene(variant, gene)) {

							BinaryInteractionQtlZscores qtlRes = binaryInteractionFile.readQtlResults(variant, gene);
							for (int j = 0; j < binaryInteractionFile.getCohortCount(); ++j) {
								sampleCountsQtl[i] = qtlRes.getSampleCounts()[j];
								zscoresQtl[i] = qtlRes.getZscores()[j];
								++i;
							}

						} else {
							for (int j = 0; j < binaryInteractionFile.getCohortCount(); ++j) {
								sampleCountsQtl[i] = -1;
								zscoresQtl[i] = Double.NaN;
								++j;
							}
						}


					}

					double metaZscore = Double.NaN; //TODO

					output.setQtlResults(variant, gene, new BinaryInteractionQtlZscores(zscoresQtl, sampleCountsQtl, metaZscore));

				}

				for (String covariate : covariates) {

						int[] sampleCountsInteraction = new int[cohortCount];

						int i = 0;
						for (BinaryInteractionFile binaryInteractionFile : binaryInteractionFiles) {

							if (binaryInteractionFile.containsInteraction(variant, gene, covariate)){

								BinaryInteractionZscores interactionRes = binaryInteractionFile.readInteractionResults(variant, gene, covariate);
								
								for (int j = 0; j < binaryInteractionFile.getCohortCount(); ++j) {
									sampleCountsInteraction[i] = interactionRes.getSamplesInteractionCohort()[j];
									++i;
								}

							} else {
								for (int j = 0; j < binaryInteractionFile.getCohortCount(); ++j) {
									sampleCountsInteraction[i] = -1;
									++j;
								}
							}


						}

						double metaZscore = Double.NaN; //TODO

						//output.setInteractionResults(variant, gene, covariate, new BinaryInteractionZscores);

					}

				}

			}

		}


	}

	private static class VariantGene {

		private final String variantName;
		private final String geneName;

		public VariantGene(String variantName, String geneName) {
			this.variantName = variantName;
			this.geneName = geneName;
		}

		public String getVariantName() {
			return variantName;
		}

		public String getGeneName() {
			return geneName;
		}

		@Override
		public int hashCode() {
			return this.variantName != null ? this.variantName.hashCode() : 0;
		}

		@Override
		public boolean equals(Object obj) {
			if (obj == null) {
				return false;
			}
			if (getClass() != obj.getClass()) {
				return false;
			}
			final VariantGene other = (VariantGene) obj;
			if ((this.variantName == null) ? (other.variantName != null) : !this.variantName.equals(other.variantName)) {
				return false;
			}
			if ((this.geneName == null) ? (other.geneName != null) : !this.geneName.equals(other.geneName)) {
				return false;
			}
			return true;
		}
	}
}