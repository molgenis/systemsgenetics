package eqtlmappingpipeline.binaryInteraction;

import eqtlmappingpipeline.Main;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.molgenis.genotype.Allele;
import umcg.genetica.io.binInteraction.BinaryInteractionCohort;
import umcg.genetica.io.binInteraction.BinaryInteractionFile;
import umcg.genetica.io.binInteraction.BinaryInteractionFileCreator;
import umcg.genetica.io.binInteraction.BinaryInteractionFileException;
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

	public static void main(String[] args) throws UnsupportedEncodingException, IOException, FileNotFoundException, BinaryInteractionFileException {

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
				} else {

					BinaryInteractionVariantCreator metaVariant = variants.get(variant.getName());

					if (!(metaVariant.getRefAllele() == variant.getRefAllele() && metaVariant.getAltAllele() == variant.getAltAllele())
							&& !(metaVariant.getRefAllele() == variant.getAltAllele() && metaVariant.getAltAllele() == variant.getRefAllele())) {
						System.err.println("Error: different alleles detected for variant: " + variant.getName());
						System.exit(1);
						return;
					}

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

		for (Map.Entry<String, BinaryInteractionVariantCreator> variantEntry : variants.entrySet()) {

			String variantName = variantEntry.getKey();
			BinaryInteractionVariantCreator variant = variantEntry.getValue();
			Allele assessedAllele = variant.getAltAllele();

			for (String gene : genes.keySet()) {

				if (allStoreQtl) {

					final int[] sampleCountsQtl = new int[cohortCount];
					final double[] zscoresQtl = new double[cohortCount];

					int i = 0;
					for (BinaryInteractionFile binaryInteractionFile : binaryInteractionFiles) {

						if (binaryInteractionFile.containsVariantGene(variantName, gene)) {

							boolean swap = binaryInteractionFile.getVariant(variantName).getAltAllele() != assessedAllele;

							BinaryInteractionQtlZscores qtlRes = binaryInteractionFile.readQtlResults(variantName, gene);
							for (int j = 0; j < binaryInteractionFile.getCohortCount(); ++j) {
								sampleCountsQtl[i] = qtlRes.getSampleCounts()[j];
								zscoresQtl[i] = qtlRes.getZscores()[j];
								if (swap) {
									zscoresQtl[i] *= -1;
								}
								++i;
							}

						} else {
							for (int j = 0; j < binaryInteractionFile.getCohortCount(); ++j) {
								sampleCountsQtl[i] = 0;
								zscoresQtl[i] = Double.NaN;
								++j;
							}
						}


					}

					double metaZscore;
					try {
						metaZscore = weightedZscore(zscoresQtl, sampleCountsQtl);
					} catch (MetaZscoreException ex) {
						if(ex.specifiedI()){
							System.err.println("Error calculating QTL meta Z-score for: " + variantName + "-" + gene + ". Problem with: " + cohorts.values().toArray(new String[cohortCount])[ex.getI()] + " error:" + ex.getMessage());
							System.exit(1);
							return;
						} else {
							System.err.println("Error calculating QTL meta Z-score for: " + variantName + "-" + gene + ". Error:" + ex.getMessage());
							System.exit(1);
							return;
						}
					}

					output.setQtlResults(variantName, gene, new BinaryInteractionQtlZscores(zscoresQtl, sampleCountsQtl, metaZscore));

				}

				for (String covariate : covariates) {

					final int[] sampleCountsInteraction = new int[cohortCount];
					final double[] zscoreSnpCohort = new double[cohortCount];
					final double[] zscoreCovariateCohort = new double[cohortCount];
					final double[] zscoreInteractionCohort = new double[cohortCount];
					final double[] rSquaredCohort = new double[cohortCount];
					final double[] zscoreInteractionFlippedCohort = new double[cohortCount];

					int i = 0;
					for (BinaryInteractionFile binaryInteractionFile : binaryInteractionFiles) {

						if (binaryInteractionFile.containsInteraction(variantName, gene, covariate)) {

							BinaryInteractionZscores interactionRes = binaryInteractionFile.readInteractionResults(variantName, gene, covariate);
							
							boolean swap = binaryInteractionFile.getVariant(variantName).getAltAllele() != assessedAllele;

							for (int j = 0; j < binaryInteractionFile.getCohortCount(); ++j) {
								sampleCountsInteraction[i] = interactionRes.getSamplesInteractionCohort()[j];
								zscoreSnpCohort[i] = interactionRes.getZscoreSnpCohort()[j];
								zscoreCovariateCohort[i] = interactionRes.getZscoreCovariateCohort()[j];
								zscoreInteractionCohort[i] = interactionRes.getZscoreInteractionCohort()[j];
								rSquaredCohort[i] = interactionRes.getrSquaredCohort()[j];
								zscoreInteractionFlippedCohort[i] = interactionRes.getZscoreInteractionFlippedCohort()[j];
								if (swap) {
									zscoreSnpCohort[i] *= -1;
									zscoreInteractionCohort[i] *= -1;
									zscoreInteractionFlippedCohort[i] *= -1;
								}
								++i;
							}

						} else {
							for (int j = 0; j < binaryInteractionFile.getCohortCount(); ++j) {
								sampleCountsInteraction[i] = 0;
								zscoreSnpCohort[i] = Double.NaN;
								zscoreCovariateCohort[i] = Double.NaN;
								zscoreInteractionCohort[i] = Double.NaN;
								rSquaredCohort[i] = Double.NaN;
								zscoreInteractionFlippedCohort[i] = Double.NaN;
								++i;
							}
						}


					}

					final double zscoreSnpMeta = Double.NaN; // TODO
					final double zscoreCovariateMeta = Double.NaN; // TODO
					final double zscoreInteractionMeta = Double.NaN; // TODO
					final double zscoreInteractionFlippedMeta;
					if (allStoreFlipped) {
						zscoreInteractionFlippedMeta = Double.NaN; // TODO
					} else {
						zscoreInteractionFlippedMeta = Double.NaN;
					}

					output.setInteractionResults(variantName, gene, covariate, new BinaryInteractionZscores(sampleCountsInteraction, zscoreSnpCohort, zscoreCovariateCohort, zscoreInteractionCohort, rSquaredCohort, zscoreInteractionFlippedCohort, zscoreSnpMeta, zscoreCovariateMeta, zscoreInteractionMeta, zscoreInteractionFlippedMeta));

				}

			}

		}

		output.close();

	}
	
	protected static double weightedZscore(double[] zscores, int[] samples) throws MetaZscoreException{
		
		double numerator = 0;
		double denominator = 0;
		
		for(int i = 0 ; i < zscores.length ; ++i){
			if(samples[i] == 0){
				continue;
			} else if(samples[i] < 0){
				throw new MetaZscoreException(i, "Sample count < 0");
			} else if(Double.isNaN(zscores[i])) {
				throw new MetaZscoreException(i, "Z-score = NaN");
			} else if(Double.isInfinite(zscores[i])){
				throw new MetaZscoreException(i, "Z-score = Inf");
			} else {
				numerator += zscores[i] * samples[i];
				denominator += samples[i] * samples[i];
			}
		}
		
		if(denominator < 1){
			throw new MetaZscoreException("No samples included");
		}
		
		return numerator / Math.sqrt(denominator);
		
	}
	
	private static class MetaZscoreException extends Exception {

		private final int i;

		public MetaZscoreException(String message) {
			super(message);
			this.i = -1;
		}
		
		public MetaZscoreException(int i, String message) {
			super(message);
			this.i = i;
		}
		
		public boolean specifiedI(){
			return i >= 0;
		}

		public int getI() {
			return i;
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