/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.transCorrelationQtl;

import au.com.bytecode.opencsv.CSVReader;
import au.com.bytecode.opencsv.CSVWriter;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.jet.random.tdouble.StudentT;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.math3.linear.SingularMatrixException;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.sampleFilter.SampleIdIncludeFilter;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variantFilter.VariantIdIncludeFilter;
import umcg.genetica.collections.intervaltree.NamedGenomicRange;
import umcg.genetica.collections.intervaltree.PerChrIntervalTree;
import umcg.genetica.containers.Pair;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class transCorrelationQtl {

	static final Pattern TAB_PATTERN = Pattern.compile("\\t");
	static final String DEFAULT_CALL_P = "0.7";
	private static final String ENCODING = "ISO-8859-1";
	private static final Options OPTIONS;
	private static final cern.jet.random.tdouble.engine.DoubleRandomEngine randomEngine = new cern.jet.random.tdouble.engine.DRand();
	private static final Pair<Double, Double> NAN_PAIR = new Pair<Double, Double>(Double.NaN, Double.NaN);
	private static final String HEADER
			= "  /---------------------------------------\\\n"
			+ "  |         Trans correlation QTL         |\n"
			+ "  |                                       |\n"
			+ "  |             Patrick Deelen            |\n"
			+ "  |        patrickdeelen@gmail.com        |\n"
			+ "  |                                       |\n"
			+ "  |     Genomics Coordination Center      |\n"
			+ "  |        Department of Genetics         |\n"
			+ "  |  University Medical Center Groningen  |\n"
			+ "  \\---------------------------------------/";

	static {

		OPTIONS = new Options();

		OptionBuilder.withArgName("basePath");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Gentoype data 1");
		OptionBuilder.withLongOpt("genotypeData");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("g"));

		OptionBuilder.withArgName("type");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("The input data 1 type.\n"
				+ "* PED_MAP - plink PED MAP files.\n"
				+ "* PLINK_BED - plink BED BIM FAM files.\n"
				+ "* VCF - bgziped vcf with tabix index file\n"
				+ "* VCFFOLDER - matches all bgziped vcf files + tabix index in a folder\n"
				+ "* SHAPEIT2 - shapeit2 phased haplotypes .haps & .sample\n"
				+ "* GEN - Oxford .gen & .sample\n"
				+ "* GENFOLDER - matches all Oxford .gen & .sample\n"
				+ "* TRITYPER - TriTyper format folder");
		OptionBuilder.withLongOpt("genotypeType");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("G"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("GTE file");
		OptionBuilder.withLongOpt("GTE");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("gte"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Expression matrix");
		OptionBuilder.withLongOpt("expression");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("e"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Trans eqtl mapping results");
		OptionBuilder.withLongOpt("transQtls");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("t"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Gene mapping file");
		OptionBuilder.withLongOpt("geneMapping");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("gm"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Output file");
		OptionBuilder.withLongOpt("output");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("o"));

	}

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException, Exception {

		System.out.println(HEADER);
		System.out.println();
		System.out.flush(); //flush to make sure header is before errors
		try {
			Thread.sleep(25); //Allows flush to complete
		} catch (InterruptedException ex) {
		}

		final String genotypeDataType;
		final String genotypeDataPath;
		final String gtePath;
		final String expressionPath;
		final String outputPath;
		final String transQtlPath;
		final String geneMappingFile;
		final int window = 100000;

		try {
			final CommandLine commandLine = new PosixParser().parse(OPTIONS, args, false);

			genotypeDataType = commandLine.getOptionValue("G");
			genotypeDataPath = commandLine.getOptionValue("g");
			gtePath = commandLine.getOptionValue("gte");
			expressionPath = commandLine.getOptionValue("e");
			outputPath = commandLine.getOptionValue("o");
			transQtlPath = commandLine.getOptionValue("t");
			geneMappingFile = commandLine.getOptionValue("gm");

		} catch (ParseException ex) {
			System.err.println("Invalid command line arguments: ");
			System.err.println(ex.getMessage());
			System.err.println();
			new HelpFormatter().printHelp(" ", OPTIONS);
			System.exit(1);
			return;
		}

		System.out.println("Genotypes: " + genotypeDataPath);
		System.out.println("Genotypes type: " + genotypeDataType);
		System.out.println("Expresison: " + expressionPath);
		System.out.println("Gte: " + gtePath);
		System.out.println("Trans eQTLs: " + transQtlPath);
		System.out.println("Gene mapping file: " + geneMappingFile);
		System.out.println("Window: " + window);
		System.out.println("Output: " + outputPath);
		System.out.println();
		System.out.println("---------------------------------");
		System.out.println();

		final HashMap<String, HashSet<Qtl>> transEqtls = loadTransEqtls(transQtlPath);
		System.out.println("Trans SNPs: " + transEqtls.size());

		final PerChrIntervalTree<NamedGenomicRange> geneMappings = loadGeneMappings(geneMappingFile, window);
		System.out.println("Genes with mapping: " + geneMappings.size());

		final Map<String, String> gteMap = readGte(gtePath);
		System.out.println("Samples in GTE: " + gteMap.size());

		final RandomAccessGenotypeData genotypeData = RandomAccessGenotypeDataReaderFormats.valueOfSmart(genotypeDataType.toUpperCase()).createFilteredGenotypeData(genotypeDataPath, 1000, new VariantIdIncludeFilter(transEqtls.keySet()), new SampleIdIncludeFilter(gteMap.keySet()));
		final String[] genotypeSamples = genotypeData.getSampleNames();
		final HashMap<String, GeneticVariant> variantMap = genotypeData.getVariantIdMap();
		System.out.println("Genotyped samples loaded: " + genotypeSamples.length);
		System.out.println("Genotyped snps loaded: " + variantMap.size());

		final DoubleMatrixDataset<String, String> expressionDataAll = DoubleMatrixDataset.loadDoubleData(expressionPath);
		System.out.println("Expression samples loaded: " + expressionDataAll.columns());

		//Put expression data in same order as genotypes
		LinkedHashSet<String> expressionSamples = new LinkedHashSet<>();
		for (String genotypeSample : genotypeSamples) {
			String expressionSample = gteMap.get(genotypeSample);
			expressionSamples.add(expressionSample);
			if (!expressionDataAll.containsCol(gteMap.get(genotypeSample))) {
				System.out.println("Could not find expression data for: " + genotypeSample + " / " + expressionSample);
			}
		}
		final DoubleMatrixDataset<String, String> expressionData = expressionDataAll.viewColSelection(expressionSamples);
		expressionData.setColObjects(Arrays.asList(genotypeData.getSampleNames()));

		final DoubleMatrixDataset<String, String> expressionDataForceNormal = expressionData.createRowForceNormalDuplicate();
		
		CSVWriter outputWriter = new CSVWriter(new FileWriter(outputPath), '\t', CSVWriter.NO_QUOTE_CHARACTER);

		final String[] outputLine = new String[8];
		int c = 0;
		outputLine[c++] = "TransSnp";
		outputLine[c++] = "Z-score";
		outputLine[c++] = "NearbyGene";
		outputLine[c++] = "TargetGene";
		outputLine[c++] = "InteractionP";
		outputLine[c++] = "InteractionZ";
		outputLine[c++] = "NumberOfSamples";
		outputLine[c++] = "RefAllele";

		StudentT tDistColt = new cern.jet.random.tdouble.StudentT(genotypeSamples.length - 4, randomEngine);

		outputWriter.writeNext(outputLine);

		for (Map.Entry<String, HashSet<Qtl>> transSnpEntry : transEqtls.entrySet()) {
			final String snp = transSnpEntry.getKey();
			final HashSet<Qtl> snpTransEffects = transSnpEntry.getValue();

			GeneticVariant variant = variantMap.get(snp);

			if (variant == null) {
				System.err.println("Skipping: " + snp);
				continue;
			}

			//List<Alleles> sampleGenotypes = variant.getSampleVariants();
			float[] sampleDosage = variant.getSampleDosages();

			if (variant.getMinorAlleleFrequency() < 0.05) {
				System.err.println("Skipping: " + snp);
				continue;
			}

			List<NamedGenomicRange> genesNearTransSnp = geneMappings.searchPosition(variant.getSequenceName(), variant.getStartPos());

			for (Qtl transEffect : snpTransEffects) {

				DoubleMatrix1D targetGeneExpression = expressionDataForceNormal.getRow(transEffect.getTrait());

				double[] targetExpresison = targetGeneExpression.toArray();
				double[][] modelParameters = new double[genotypeSamples.length][];

				for (NamedGenomicRange geneNearTransSnp : genesNearTransSnp) {

					if (!expressionDataForceNormal.containsRow(geneNearTransSnp.getName())) {
						continue;
					}

					DoubleMatrix1D nearbyGeneExpression = expressionDataForceNormal.getRow(geneNearTransSnp.getName());

					for (int s = 0; s < genotypeSamples.length; ++s) {

						double[] x = new double[3];
						x[0] = nearbyGeneExpression.getQuick(s);
						x[1] = sampleDosage[s];
						x[2] = nearbyGeneExpression.getQuick(s) * sampleDosage[s];

						modelParameters[s] = x;

					}

					double pValueInteraction;
					double zScoreInteraction;

					try {
						OLSMultipleLinearRegression regression = new OLSMultipleLinearRegression();
						regression.newSampleData(targetExpresison, modelParameters);
						double[] betas = regression.estimateRegressionParameters();
						double[] betaSes = regression.estimateRegressionParametersStandardErrors();

						Pair<Double, Double> pair = convertBetaToP(betas[3], betaSes[3], tDistColt);
						pValueInteraction = pair.getLeft();
						zScoreInteraction = pair.getRight();
					} catch (SingularMatrixException ex) {
						pValueInteraction = Double.NaN;
						zScoreInteraction = Double.NaN;
					}

					c = 0;
					outputLine[c++] = snp;
					outputLine[c++] = String.valueOf(transEffect.getZscore());
					outputLine[c++] = geneNearTransSnp.getName();
					outputLine[c++] = transEffect.getTrait();
					outputLine[c++] = String.valueOf(pValueInteraction);
					outputLine[c++] = String.valueOf(zScoreInteraction);
					outputLine[c++] = String.valueOf(genotypeSamples.length);
					outputLine[c++] = variant.getRefAllele().getAlleleAsString();
					outputWriter.writeNext(outputLine);

				}
			}

		}

		outputWriter.close();

	}

	public static PerChrIntervalTree<NamedGenomicRange> loadGeneMappings(final String geneMappingFile, final int window) throws IOException, NumberFormatException, Exception {
		CSVReader geneMapReader = new CSVReader(new InputStreamReader(new FileInputStream(geneMappingFile), ENCODING), '\t', '\0', 1);
		String[] nextLine;
		HashMap<String, ArrayList<NamedGenomicRange>> genes = new HashMap<>();
		while ((nextLine = geneMapReader.readNext()) != null) {
			String name = nextLine[0];
			String chr = nextLine[1];
			int start = Integer.valueOf(nextLine[2]);
			int stop = Integer.valueOf(nextLine[3]);

			ArrayList<NamedGenomicRange> chrGenes = genes.getOrDefault(chr, new ArrayList<>());
			chrGenes.add(new NamedGenomicRange(name, chr, start - window < 0 ? 0 : start - window, stop + window));
			genes.putIfAbsent(chr, chrGenes);

		}
		PerChrIntervalTree<NamedGenomicRange> geneMappings = new PerChrIntervalTree(NamedGenomicRange.class);
		for (Map.Entry<String, ArrayList<NamedGenomicRange>> entry : genes.entrySet()) {
			geneMappings.addChrElements(entry.getKey(), entry.getValue());
		}
		return geneMappings;
	}

//	public static HashMap<String, HashSet<EQTL>> loadTransEqtls(final String transQtlPath) throws IOException {
//		final QTLTextFile transEqtlFile = new QTLTextFile(transQtlPath, false);
//		final HashMap<String, HashSet<EQTL>> transEqtls = new HashMap<>();
//		Iterator<EQTL> transEqtlFileIterator = transEqtlFile.getEQtlIterator();
//		while (transEqtlFileIterator.hasNext()) {
//			EQTL eqtl = transEqtlFileIterator.next();
//			HashSet<EQTL> snpTransEffects = transEqtls.getOrDefault(eqtl.getRsName(), new HashSet());
//			snpTransEffects.add(eqtl);
//			transEqtls.putIfAbsent(eqtl.getRsName(), snpTransEffects);
//		}
//		return transEqtls;
//	}
	public static HashMap<String, HashSet<Qtl>> loadTransEqtls(final String transQtlPath) throws IOException {

		final HashMap<String, HashSet<Qtl>> transEqtls = new HashMap<>();

		CSVReader transFileReader = new CSVReader(new InputStreamReader(new FileInputStream(transQtlPath), ENCODING), '\t', '\0', 0);
		String[] nextLine;
		while ((nextLine = transFileReader.readNext()) != null) {

			Qtl eqtl = new Qtl(nextLine[3], nextLine[4], Double.valueOf(nextLine[5]));

			HashSet<Qtl> snpTransEffects = transEqtls.getOrDefault(eqtl.getSnp(), new HashSet());
			snpTransEffects.add(eqtl);
			transEqtls.putIfAbsent(eqtl.getSnp(), snpTransEffects);

		}

		return transEqtls;
	}

	public static Map<String, String> readGte(final String gtePath) throws IOException {
		final HashMap<String, String> genotypeToExpression = new HashMap<>();
		CSVReader gteFileReader = new CSVReader(new InputStreamReader(new FileInputStream(gtePath), ENCODING), '\t', '\0', 0);
		String[] nextLine;
		while ((nextLine = gteFileReader.readNext()) != null) {

			final String genotypeId = nextLine[0];
			final String expressionId = nextLine[1];

			genotypeToExpression.put(genotypeId, expressionId);

		}
		return Collections.unmodifiableMap(genotypeToExpression);
	}

	private static Pair<Double, Double> convertBetaToP(double beta, double se, StudentT tDistColt) {

		if (Double.isNaN(beta)) {
			return NAN_PAIR;
		}

		double t = beta / se;
		double p = 1;
		double z = 0;
		if (t < 0) {
			p = tDistColt.cdf(t);
			if (p < 2.0E-323) {
				p = 2.0E-323;

			}
			z = cern.jet.stat.Probability.normalInverse(p);
		} else {
			p = tDistColt.cdf(-t);
			if (p < 2.0E-323) {
				p = 2.0E-323;

			}
			z = -cern.jet.stat.Probability.normalInverse(p);
		}
		return new Pair<Double, Double>(p, z);
	}

}
