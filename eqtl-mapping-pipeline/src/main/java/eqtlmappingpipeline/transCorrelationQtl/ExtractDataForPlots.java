/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.transCorrelatieQtl;

import au.com.bytecode.opencsv.CSVReader;
import au.com.bytecode.opencsv.CSVWriter;
import static eqtlmappingpipeline.transCorrelatieQtl.transCorrelationQtl.readGte;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.regex.Pattern;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.sampleFilter.SampleIdIncludeFilter;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variantFilter.VariantIdIncludeFilter;
import umcg.genetica.containers.Pair;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class ExtractDataForPlots {

	static final Pattern TAB_PATTERN = Pattern.compile("\\t");
	static final String DEFAULT_CALL_P = "0.7";
	private static final String ENCODING = "ISO-8859-1";
	private static final Options OPTIONS;
	private static final cern.jet.random.tdouble.engine.DoubleRandomEngine randomEngine = new cern.jet.random.tdouble.engine.DRand();
	private static final Pair<Double, Double> NAN_PAIR = new Pair<Double, Double>(Double.NaN, Double.NaN);
	private static final String HEADER
			= "  /---------------------------------------\\\n"
			+ "  |      Meta Trans correlation QTL       |\n"
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

		OptionBuilder.withArgName("strings");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("cohort names");
		OptionBuilder.withLongOpt("cohorts");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("c"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("GenotypeFolder");
		OptionBuilder.withLongOpt("g");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("g"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Output file");
		OptionBuilder.withLongOpt("output");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("o"));

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
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Interactions to extract");
		OptionBuilder.withLongOpt("interactions");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("i"));

	}

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws Exception {

		System.out.println(HEADER);
		System.out.println();
		System.out.flush(); //flush to make sure header is before errors
		try {
			Thread.sleep(25); //Allows flush to complete
		} catch (InterruptedException ex) {
		}

		final String[] cohorts;
		final String genotypeFolder;
		final String outputPath;
		final String gtePath;
		final String expressionPath;
		final String interactionsPath;

		try {
			final CommandLine commandLine = new PosixParser().parse(OPTIONS, args, false);

			cohorts = commandLine.getOptionValues("c");
			genotypeFolder = commandLine.getOptionValue("g");
			outputPath = commandLine.getOptionValue("o");
			gtePath = commandLine.getOptionValue("gte");
			expressionPath = commandLine.getOptionValue("e");
			interactionsPath = commandLine.getOptionValue("i");

		} catch (ParseException ex) {
			System.err.println("Invalid command line arguments: ");
			System.err.println(ex.getMessage());
			System.err.println();
			new HelpFormatter().printHelp(" ", OPTIONS);
			System.exit(1);
			return;
		}

		HashSet<String> snpsToExtract = new HashSet<>();
		LinkedHashSet<String> genesToExtract = new LinkedHashSet<>();
		
		CSVReader geneMapReader = new CSVReader(new InputStreamReader(new FileInputStream(interactionsPath), ENCODING), '\t', '\0', 1);
		String[] nextLine;
		while ((nextLine = geneMapReader.readNext()) != null) {
			String snp = nextLine[0];
			String localGene = nextLine[2];
			String targetGene = nextLine[3];
			
			snpsToExtract.add(snp);
			genesToExtract.add(localGene);
			genesToExtract.add(targetGene);
			
		}
		
		final Map<String, String> gteMap = readGte(gtePath);
		System.out.println("Samples in GTE: " + gteMap.size());

		HashMap<String, RandomAccessGenotypeData> genotypeDatasets = new HashMap<>();
		HashMap<String, HashMap<String, GeneticVariant>> perCohortVariantMap = new HashMap<>();
		ArrayList<String> genotypedSamples = new ArrayList<>();
		
		
		for(String cohort : cohorts){
			final RandomAccessGenotypeData genotypeData = RandomAccessGenotypeDataReaderFormats.TRITYPER.createFilteredGenotypeData(genotypeFolder + "/" + cohort + "/combined/", 1000, new VariantIdIncludeFilter(snpsToExtract), new SampleIdIncludeFilter(gteMap.keySet()));
			perCohortVariantMap.put(cohort, genotypeData.getVariantIdMap());
			genotypeDatasets.put(cohort, genotypeData);
			for(String sample : genotypeData.getSampleNames()){
				genotypedSamples.add(sample);
			}
			
		}
		
		System.out.println("Total number of genotyped samples: " + genotypedSamples.size());
		
		final DoubleMatrixDataset<String, String> expressionDataAll = DoubleMatrixDataset.loadDoubleData(expressionPath);
		System.out.println("Expression samples loaded: " + expressionDataAll.columns());

		//Put expression data in same order as genotypes
		LinkedHashSet<String> expressionSamples = new LinkedHashSet<>();
		for (String genotypeSample : genotypedSamples) {
			String expressionSample = gteMap.get(genotypeSample);
			expressionSamples.add(expressionSample);
			if (!expressionDataAll.containsCol(gteMap.get(genotypeSample))) {
				System.out.println("Could not find expression data for: " + genotypeSample + " / " + expressionSample);
			}
		}
		final DoubleMatrixDataset<String, String> expressionData = expressionDataAll.viewColSelection(expressionSamples).viewRowSelection(genesToExtract);
		
		expressionData.setColObjects(genotypedSamples);

		//final DoubleMatrixDataset<String, String> expressionDataForceNormal = expressionData.createRowForceNormalDuplicate();
		
		final DoubleMatrixDataset<String, String> expressionDataForceNormal = expressionData;
		
		expressionDataForceNormal.save(outputPath + "_expression.txt");
		
		CSVWriter outputWriter = new CSVWriter(new FileWriter(outputPath + "_genotypes.txt"), '\t', CSVWriter.NO_QUOTE_CHARACTER);

		final String[] outputLine = new String[genotypedSamples.size()+1];
		outputLine[0] = "SNP";
		for(int s = 1; s < genotypedSamples.size() ; ++s){
			outputLine[s] = genotypedSamples.get(s);
		}
		outputWriter.writeNext(outputLine);
		
		for(String snp : snpsToExtract){
			int s = 0;
			outputLine[s++] = snp;
			for(String cohort : cohorts){
				GeneticVariant variant = perCohortVariantMap.get(cohort).get(snp);
				int numberOfSamples = genotypeDatasets.get(cohort).getSampleNames().length;
				
				if(variant == null || variant.getMinorAlleleFrequency() < 0.05){
					
					for(int i = 0 ; i < numberOfSamples ; ++i){
						outputLine[s++] = "NA";
					}
					
				} else {
					
					for(Alleles a : variant.getSampleVariants()){
						outputLine[s++] = a.toString();
					}
					
				}
				
				
			}
			outputWriter.writeNext(outputLine);
		}
		outputWriter.close();
	}

}
