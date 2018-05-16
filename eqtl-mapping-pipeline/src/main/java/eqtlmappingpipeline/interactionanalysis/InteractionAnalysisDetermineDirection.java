package eqtlmappingpipeline.interactionanalysis;

import au.com.bytecode.opencsv.CSVReader;
import au.com.bytecode.opencsv.CSVWriter;
import gnu.trove.list.array.TDoubleArrayList;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.Map;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.bidimap.DualHashBidiMap;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.apache.commons.math3.stat.ranking.NaNStrategy;
import org.apache.commons.math3.stat.ranking.NaturalRanking;
import org.apache.commons.math3.stat.ranking.RankingAlgorithm;
import org.apache.log4j.Logger;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.GenotypeInfo;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.multipart.IncompatibleMultiPartGenotypeDataException;
import org.molgenis.genotype.tabix.TabixFileNotFoundException;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author Patrick Deelen
 */
public class InteractionAnalysisDetermineDirection {

	private final RandomAccessGenotypeData genotypeData;
	private final DoubleMatrixDataset<String, String> expressionData;
	private final DoubleMatrixDataset<String, String> covariatesData;
	private final BidiMap<String, String> gte;
	private final HashMap<String, GeneticVariant> variantIdMap;
	private static final RankingAlgorithm COV_RANKER = new NaturalRanking(NaNStrategy.FAILED, new Well19937c(1));
	private static final SpearmansCorrelation spearmanCalculator = new SpearmansCorrelation();
	private static final Options OPTIONS;
	private static Logger LOGGER;

	static {

		LOGGER = Logger.getLogger(GenotypeInfo.class);

		OPTIONS = new Options();

		OptionBuilder.withArgName("basePath");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("The genotype");
		OptionBuilder.withLongOpt("genotypes");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("g"));

		OptionBuilder.withArgName("format");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("The genotype data format. If not defined will attempt to automatically select the first matching dataset on the specified path\n"
				+ "* PED_MAP - plink PED MAP files.\n"
				+ "* PLINK_BED - plink BED BIM FAM files.\n"
				+ "* VCF - bgziped vcf with tabix index file\n"
				+ "* VCFFOLDER - matches all bgziped vcf files + tabix index in a folder\n"
				+ "* SHAPEIT2 - shapeit2 phased haplotypes .haps & .sample\n"
				+ "* GEN - Oxford .gen & .sample\n"
				+ "* TRITYPER - TriTyper format folder");
		OptionBuilder.withLongOpt("genotypesFormat");
		OPTIONS.addOption(OptionBuilder.create("G"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Expression data");
		OptionBuilder.withLongOpt("expression");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("e"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Covariate data");
		OptionBuilder.withLongOpt("covariates");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("c"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Genotype to expression coupling");
		OptionBuilder.withLongOpt("gte");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("gte"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Query variant <tab> gene <tab> covariate <tab> assessedAllele. No header");
		OptionBuilder.withLongOpt("query");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("q"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Output file");
		OptionBuilder.withLongOpt("output");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("o"));

		OptionBuilder.withArgName("double");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Fraction of tail of either end of covarate to use.");
		OptionBuilder.withLongOpt("fraction");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("f"));

	}

	public static void main(String[] args) throws IOException {

		CommandLineParser parser = new PosixParser();
		final CommandLine commandLine;
		try {
			commandLine = parser.parse(OPTIONS, args, false);
		} catch (ParseException ex) {
			System.err.println("Invalid command line arguments: " + ex.getMessage());
			System.err.println();
			new HelpFormatter().printHelp(" ", OPTIONS);
			System.exit(1);
			return;
		}

		final String[] genotypePath = commandLine.getOptionValues("g");
		final RandomAccessGenotypeDataReaderFormats genotypeFormat;

		try {
			if (commandLine.hasOption("G")) {
				genotypeFormat = RandomAccessGenotypeDataReaderFormats.valueOf(commandLine.getOptionValue("G").toUpperCase());
			} else {
				if (genotypePath[0].endsWith(".vcf")) {
					System.err.println("Only vcf.gz is supported. Please see manual on how to do create a vcf.gz file.");
					System.exit(1);
					return;
				}
				try {
					genotypeFormat = RandomAccessGenotypeDataReaderFormats.matchFormatToPath(genotypePath[0]);
				} catch (GenotypeDataException e) {
					System.err.println("Unable to determine input 1 type based on specified path. Please specify --G");
					System.exit(1);
					return;
				}
			}
		} catch (IllegalArgumentException e) {
			System.err.println("Error parsing --G \"" + commandLine.getOptionValue("G") + "\" is not a valid input data format");
			System.exit(1);
			return;
		}

		final String expressionDataPath = commandLine.getOptionValue("e");
		final String covariateDataPath = commandLine.getOptionValue("c");
		final String gtePath = commandLine.getOptionValue("gte");
		final String queryPath = commandLine.getOptionValue("q");
		final String outputPath = commandLine.getOptionValue("o");
		final double fractionToUse = Double.parseDouble(commandLine.getOptionValue("f"));

		System.out.println("Genotype data: " + genotypePath);
		System.out.println("Genotype data format: " + genotypeFormat);
		System.out.println("Expression data: " + expressionDataPath);
		System.out.println("Covariate data: " + covariateDataPath);
		System.out.println("Gte data: " + gtePath);
		System.out.println("Query: " + queryPath);
		System.out.println("Output: " + outputPath);
		System.out.println("Outer fractions to use: " + fractionToUse);

		final RandomAccessGenotypeData genotypeData;

		try {
			genotypeData = genotypeFormat.createFilteredGenotypeData(genotypePath, 100, null, null, null, 0.8);
		} catch (TabixFileNotFoundException e) {
			LOGGER.fatal("Tabix file not found for input data at: " + e.getPath() + "\n"
					+ "Please see README on how to create a tabix file");
			System.exit(1);
			return;
		} catch (IOException e) {
			LOGGER.fatal("Error reading input data: " + e.getMessage(), e);
			System.exit(1);
			return;
		} catch (IncompatibleMultiPartGenotypeDataException e) {
			LOGGER.fatal("Error combining the impute genotype data files: " + e.getMessage(), e);
			System.exit(1);
			return;
		} catch (GenotypeDataException e) {
			LOGGER.fatal("Error reading input data: " + e.getMessage(), e);
			System.exit(1);
			return;
		}

		System.out.println("Genotype data loaded for " + genotypeData.getSampleNames().length + " individuals");

		final DoubleMatrixDataset<String, String> expressionData = DoubleMatrixDataset.loadDoubleTextData(expressionDataPath, '\t');

		System.out.println("Loaded expression data for: " + expressionData.rows() + " genes and " + expressionData.columns() + " individuals");

		final DoubleMatrixDataset<String, String> covariatesData = DoubleMatrixDataset.loadDoubleTextData(covariateDataPath, '\t');

		System.out.println("Loaded covariate data for: " + expressionData.rows() + " genes and " + expressionData.columns() + " individuals");

		final BidiMap<String, String> gte = loadGte(gtePath);

		InteractionAnalysisDetermineDirection directionTool = new InteractionAnalysisDetermineDirection(genotypeData, expressionData, covariatesData, gte);

		CSVReader reader = new CSVReader(new FileReader(queryPath), '\t', '\0', 1);
		CSVWriter writer = new CSVWriter(new FileWriter(outputPath), '\t', CSVWriter.NO_QUOTE_CHARACTER);

		String[] outputLine = new String[6];
		int c = 0;
		outputLine[c++] = "variant";
		outputLine[c++] = "gene";
		outputLine[c++] = "covariate";
		outputLine[c++] = "assessedAllele";
		outputLine[c++] = "rhoLow";
		outputLine[c++] = "rhoHigh";
		writer.writeNext(outputLine);

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {
			
			final String variant = nextLine[0];
			final String gene = nextLine[1];
			final String covariate = nextLine[2];
			final Allele assessedAllele = Allele.create(nextLine[3]);

			final EffectDiffResult effectDiff = directionTool.calculateEffectDifference(variant, gene, covariate, assessedAllele, fractionToUse);

			c = 0;
			outputLine[c++] = variant;
			outputLine[c++] = gene;
			outputLine[c++] = covariate;
			outputLine[c++] = assessedAllele.getAlleleAsString();
			outputLine[c++] = String.valueOf(effectDiff.getRhoLow());
			outputLine[c++] = String.valueOf(effectDiff.getRhoHigh());
			writer.writeNext(outputLine);

		}
		writer.close();
		reader.close();
		System.out.println("Done");

	}

	private static BidiMap<String, String> loadGte(String gtePath) throws IOException {

		BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(gtePath), "UTF-8"));

		String line;

		BidiMap<String, String> gte = new DualHashBidiMap<String, String>();

		while ((line = reader.readLine()) != null) {
			String[] elements = StringUtils.split(line, '\t');
			if (elements.length != 2) {
				throw new RuntimeException("Error in GTE file line: " + line);
			}
			gte.put(elements[0], elements[1]);
		}

		return gte;
	}

	public InteractionAnalysisDetermineDirection(RandomAccessGenotypeData genotypeData, DoubleMatrixDataset<String, String> expressionData, DoubleMatrixDataset<String, String> covariatesData, BidiMap<String, String> gte) {
		this.genotypeData = genotypeData;
		this.expressionData = expressionData;
		this.covariatesData = covariatesData;
		this.gte = gte;
		this.variantIdMap = genotypeData.getVariantIdMap();

		HashSet<String> genotypedSamples = new HashSet<String>();
		Collections.addAll(genotypedSamples, genotypeData.getSampleNames());

		for (Iterator<Map.Entry<String, String>> it = gte.entrySet().iterator(); it.hasNext();) {
			Map.Entry<String, String> gteEntry = it.next();

			if (!genotypedSamples.contains(gteEntry.getKey())) {
				it.remove();
			}

			if (!expressionData.containsCol(gteEntry.getValue())) {
				it.remove();
			}

			if (!covariatesData.containsCol(gteEntry.getValue())) {
				it.remove();
			}

		}

		System.out.println("Samples with: genotypes, expression & covariate data: " + gte.size());

	}

	public EffectDiffResult calculateEffectDifference(String snpId, String geneName, String covariateName, Allele assessedAllele, double fractionOfSamplesPerGroup) {

		if (!variantIdMap.containsKey(snpId)) {
			return new EffectDiffResult(Double.NaN, Double.NaN);
		}

		if (!expressionData.containsRow(geneName)) {
			return new EffectDiffResult(Double.NaN, Double.NaN);
		}

		if (!covariatesData.containsRow(covariateName)) {
			return new EffectDiffResult(Double.NaN, Double.NaN);
		}

		if (fractionOfSamplesPerGroup <= 0 || fractionOfSamplesPerGroup >= 1) {
			throw new RuntimeException("Fraction must be between 0 and 1");
		}

		GeneticVariant variant = variantIdMap.get(snpId);
		Alleles variantAlleles = variant.getVariantAlleles();

		if (!variantAlleles.contains(assessedAllele)) {
			return new EffectDiffResult(Double.NaN, Double.NaN);
		}

		if (variantAlleles.getAlleleCount() != 2) {
			return new EffectDiffResult(Double.NaN, Double.NaN);
		}

		float[] dosagesAll = variant.getSampleDosages();
		String[] genotypedSamples = genotypeData.getSampleNames();

		LinkedHashSet<String> includedGenotypedSamples = new LinkedHashSet<>();
		TDoubleArrayList dosages = new TDoubleArrayList(dosagesAll.length);

		for (int i = 0; i < dosagesAll.length; ++i) {
			if (dosagesAll[i] >= 0 && gte.containsKey(genotypedSamples[i])) {
				includedGenotypedSamples.add(genotypedSamples[i]);
				dosages.add(dosagesAll[i]);
			}
		}

		System.out.println("Included samples: " + includedGenotypedSamples.size());

		double[] expressionLevels = new double[includedGenotypedSamples.size()];
		double[] covariateLevels = new double[includedGenotypedSamples.size()];

		int s = 0;
		for (String genotypeSample : includedGenotypedSamples) {
			expressionLevels[s] = expressionData.getElement(geneName, gte.get(genotypeSample));
			covariateLevels[s] = covariatesData.getElement(covariateName, gte.get(genotypeSample));
			++s;
		}

		if (assessedAllele != variantAlleles.get(0)) {
			for (int i = 0; i < dosages.size(); ++i) {
				dosages.setQuick(i, dosages.getQuick(i) * -1);
			}
		}

		double[] covariateRanks = COV_RANKER.rank(covariateLevels);

		int samplesPerGroup = (int) Math.floor(covariateRanks.length * fractionOfSamplesPerGroup);

		System.out.println("Samples per group: " + samplesPerGroup);

		double[] dosagesLow = new double[samplesPerGroup];
		double[] expressionLow = new double[samplesPerGroup];

		double[] dosagesHigh = new double[samplesPerGroup];
		double[] expressionHigh = new double[samplesPerGroup];

		for (int i = 0; i < samplesPerGroup; ++i) {
			dosagesLow[i] = dosages.get((int) covariateRanks[i]);
			expressionLow[i] = expressionLevels[(int) covariateRanks[i]];
			dosagesHigh[i] = dosages.get((int) covariateRanks[covariateRanks.length - 1 - i]);
			expressionHigh[i] = expressionLevels[(int) covariateRanks[covariateRanks.length - 1 - i]];
		}

		double rhoLow = spearmanCalculator.correlation(dosagesLow, expressionLow);
		double rhoHigh = spearmanCalculator.correlation(dosagesHigh, expressionHigh);

		System.out.println("rho low:" + rhoLow);
		System.out.println("rho high:" + rhoHigh);

		return new EffectDiffResult(rhoLow, rhoHigh);

	}

	static class EffectDiffResult {

		private final double rhoLow;
		private final double rhoHigh;

		public EffectDiffResult(double rhoLow, double rhoHigh) {
			this.rhoLow = rhoLow;
			this.rhoHigh = rhoHigh;
		}

		public double getRhoLow() {
			return rhoLow;
		}

		public double getRhoHigh() {
			return rhoHigh;
		}
	}
}
