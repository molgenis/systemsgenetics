package eqtlmappingpipeline.binaryInteraction;

import au.com.bytecode.opencsv.CSVReader;
import eqtlmappingpipeline.Main;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.UnsupportedEncodingException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.lang3.StringUtils;
import org.molgenis.genotype.Allele;
import umcg.genetica.io.binInteraction.BinaryInteractionCohort;
import umcg.genetica.io.binInteraction.BinaryInteractionFile;
import umcg.genetica.io.binInteraction.BinaryInteractionFileCreator;
import umcg.genetica.io.binInteraction.BinaryInteractionQtlZscores;
import umcg.genetica.io.binInteraction.BinaryInteractionZscores;
import umcg.genetica.io.binInteraction.gene.BinaryInteractionGeneCreator;
import umcg.genetica.io.binInteraction.variant.BinaryInteractionVariantCreator;

/**
 *
 * @author Patrick Deelen
 */
public class ConvertTextOutputToBinary {

	private static final String VERSION = Main.VERSION;
	private static final String HEADER =
			"  /---------------------------------------\\\n"
			+ "  |     Convert text output to binary     |\n"
			+ "  |                                       |\n"
			+ "  |             Patrick Deelen            |\n"
			+ "  |        patrickdeelen@gmail.com        |\n"
			+ "  |                                       |\n"
			+ "  |   Dasha Zhernakova, Marc Jan Bonder   |\n"
			+ "  |      Lude Franke, Morris Swertz       |\n"
			+ "  |                                       |\n"
			+ "  |     Genomics Coordination Center      |\n"
			+ "  |        Department of Genetics         |\n"
			+ "  |  University Medical Center Groningen  |\n"
			+ "  \\---------------------------------------/";
	private static final DateFormat DATE_TIME_FORMAT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
	private static final Date currentDataTime = new Date();
	private static final Options OPTIONS;
	private static final String ENCODING = "ISO-8859-1";
	private static final Pattern SLASH_PATTERN = Pattern.compile("/");

	static {

		OPTIONS = new Options();

		Option option;

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Input interaction file");
		OptionBuilder.withLongOpt("interaction");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("i"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Input SNPs file");
		OptionBuilder.withLongOpt("snps");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("s"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Output file");
		OptionBuilder.withLongOpt("output");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("o"));

		OptionBuilder.withArgName("string");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Cohort name");
		OptionBuilder.withLongOpt("cohort");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("c"));

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

		final File inputInteractionFile;
		final File inputSnpFile;
		final File outputFile;
		final String cohortName;

		try {
			final CommandLine commandLine = new PosixParser().parse(OPTIONS, args, false);

			inputInteractionFile = new File(commandLine.getOptionValue("i"));
			inputSnpFile = new File(commandLine.getOptionValue("s"));
			outputFile = new File(commandLine.getOptionValue("o"));
			cohortName = commandLine.getOptionValue("c");

		} catch (ParseException ex) {
			System.err.println("Invalid command line arguments: ");
			System.err.println(ex.getMessage());
			System.err.println();
			new HelpFormatter().printHelp(" ", OPTIONS);
			System.exit(1);
			return;
		}

		System.out.println("Input interaction file: " + inputInteractionFile.getAbsolutePath());
		System.out.println("Input SNP file: " + inputSnpFile.getAbsolutePath());
		System.out.println("Output file: " + outputFile.getAbsolutePath());
		System.out.println("Cohort name: " + cohortName);

		if (!inputInteractionFile.canRead()) {
			System.err.println("Cannot read input file.");
			System.exit(1);
			return;
		}

		if (!inputSnpFile.canRead()) {
			System.err.println("Cannot read input file.");
			System.exit(1);
			return;
		}

		final BufferedReader inputInteractionReaderRun1;
		if (inputInteractionFile.getName().endsWith(".gz")) {
			inputInteractionReaderRun1 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inputInteractionFile)), ENCODING));
		} else {
			inputInteractionReaderRun1 = new BufferedReader(new InputStreamReader(new FileInputStream(inputInteractionFile), ENCODING));
		}

		final BufferedReader inputInteractionReaderRun2;
		if (inputInteractionFile.getName().endsWith(".gz")) {
			inputInteractionReaderRun2 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inputInteractionFile)), ENCODING));
		} else {
			inputInteractionReaderRun2 = new BufferedReader(new InputStreamReader(new FileInputStream(inputInteractionFile), ENCODING));
		}

		final CSVReader inputSnpReader;
		if (inputSnpFile.getName().endsWith(".gz")) {
			inputSnpReader = new CSVReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inputSnpFile)), ENCODING), '\t', '\0', 1);
		} else {
			inputSnpReader = new CSVReader(new InputStreamReader(new FileInputStream(inputSnpFile), ENCODING), '\t', '\0', 1);
		}

		ArrayList<BinaryInteractionVariantCreator> variants = new ArrayList<BinaryInteractionVariantCreator>();

		String[] nextLine;
		while ((nextLine = inputSnpReader.readNext()) != null) {

			final String snpName = nextLine[0];
			final String chr = nextLine[1];
			final int pos = Integer.parseInt(nextLine[2]);
			final String[] alleles = SLASH_PATTERN.split(nextLine[3]);
			final Allele ma = Allele.create(nextLine[4]);

			Allele a1 = Allele.create(alleles[0]);
			Allele a2 = Allele.create(alleles[1]);

			if (ma == a1) {
				a1 = a2;
				a2 = ma;
			} else if (ma != a2) {
				throw new Exception("Minor allele not one of the alles for SNP: " + snpName);
			}

			variants.add(new BinaryInteractionVariantCreator(snpName, chr, pos, a1, a2));

		}
		inputSnpReader.close();

		System.out.println("Completed reading " + variants.size() + " variants");

		LinkedHashMap<String, BinaryInteractionGeneCreator> genes = new LinkedHashMap<String, BinaryInteractionGeneCreator>();
		LinkedHashSet<String> covariates = new LinkedHashSet<String>();
		LinkedHashSet<VariantGene> variantGenes = new LinkedHashSet<VariantGene>();

		String line;
		inputInteractionReaderRun1.readLine();
		int lastReport = 0;
		while ((line = inputInteractionReaderRun1.readLine()) != null) {

			nextLine = StringUtils.split(line, '\t');

			final String variantName = nextLine[0];
			final String geneName = nextLine[1];
			final String covariateName = nextLine[2];

			if (!genes.containsKey(geneName)) {
				genes.put(geneName, new BinaryInteractionGeneCreator(geneName));
			}

			covariates.add(covariateName);

			variantGenes.add(new VariantGene(variantName, geneName));

//			if (variantGenes.size() != lastReport && variantGenes.size() % 1000 == 0) {
//				lastReport = variantGenes.size();
//				System.out.println("Running first round, so far detected:");
//				System.out.println(" - " + genes.size() + " genes");
//				System.out.println(" - " + covariates.size() + " covariates");
//				System.out.println(" - " + variantGenes.size() + " variant-gene combinations");
//			}

		}

		inputInteractionReaderRun1.close();

		System.out.println("Completed first round parsing interaction text file");
		System.out.println("Detected:");
		System.out.println(" - " + genes.size() + " genes");
		System.out.println(" - " + covariates.size() + " covariates");
		System.out.println(" - " + variantGenes.size() + " variant-gene combinations");
		System.out.println();

		BinaryInteractionCohort[] cohorts = {new BinaryInteractionCohort(cohortName, -1)};

		BinaryInteractionFileCreator binaryFileCreator = new BinaryInteractionFileCreator(outputFile, variants.toArray(new BinaryInteractionVariantCreator[variants.size()]), genes.values().toArray(new BinaryInteractionGeneCreator[genes.size()]), cohorts, covariates.toArray(new String[covariates.size()]), true, false, true, true);

		for (VariantGene variantGene : variantGenes) {
			binaryFileCreator.addTestedVariantGene(variantGene.getVariantName(), variantGene.getGeneName());
		}

		binaryFileCreator.setDescription("File converted from txt result: " + inputInteractionFile.getAbsolutePath() + ". Using eQTL pipeline version: " + VERSION);

		BinaryInteractionFile binaryInteractionFile = binaryFileCreator.create();

		String lastGene = "";
		String lastVariant = "";

		inputInteractionReaderRun2.readLine();
		while ((line = inputInteractionReaderRun2.readLine()) != null) {

			nextLine = StringUtils.split(line, '\t');

			final String variantName = nextLine[0];
			final String geneName = nextLine[1];
			final String covariateName = nextLine[2];

			final double zscoreSnpCohort = Double.parseDouble(nextLine[3]);
			final double zscoreCovariateCohort = Double.parseDouble(nextLine[4]);
			final double zscoreInteractionCohort = Double.parseDouble(nextLine[5]);
			final double zscoreQtl = Double.parseDouble(nextLine[6]);
			final double zscoreInteractionFlippedCohort = Double.parseDouble(nextLine[7]);
			final int samplesInteractionCohort = Integer.parseInt(nextLine[8]);
			final double rSquaredCohort = Double.parseDouble(nextLine[9]);

			BinaryInteractionZscores zscoresInteractions = new BinaryInteractionZscores(samplesInteractionCohort, zscoreSnpCohort, zscoreCovariateCohort, zscoreInteractionCohort, rSquaredCohort, zscoreInteractionFlippedCohort);

			binaryInteractionFile.setInteractionResults(variantName, geneName, covariateName, zscoresInteractions);

			if (!lastGene.equals(geneName) || !lastVariant.equals(variantName)) {
				BinaryInteractionQtlZscores zscoresQtl = new BinaryInteractionQtlZscores(zscoreQtl, samplesInteractionCohort);
				binaryInteractionFile.setQtlResults(variantName, geneName, zscoresQtl);
				lastGene = geneName;
				lastVariant = variantName;
			}


		}

		inputInteractionReaderRun2.close();

		binaryInteractionFile.finalizeWriting();
		binaryInteractionFile.close();

		System.out.println("Completed binary file");
		System.out.println("Writen:");
		System.out.println(" - " + binaryInteractionFile.getInteractionZscoresSet() + " interaction with z-scores written. Using " + binaryInteractionFile.getInteractionWriteBufferFlushed() + " buffer flushes");
		System.out.println(" - " + binaryInteractionFile.getQtlZscoresSet() + " qtls with z-scores written. Using " + binaryInteractionFile.getQtlWriteBufferFlushed() + " buffer flushes");
		System.out.println("");
		System.out.println("Current date and time: " + DATE_TIME_FORMAT.format(new Date()));


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
			int hash = 7;
			hash = 31 * hash + (this.variantName != null ? this.variantName.hashCode() : 0);
			hash = 31 * hash + (this.geneName != null ? this.geneName.hashCode() : 0);
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
