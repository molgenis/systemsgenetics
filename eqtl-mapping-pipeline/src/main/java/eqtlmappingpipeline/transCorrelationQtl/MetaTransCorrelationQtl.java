/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.transCorrelationQtl;

import au.com.bytecode.opencsv.CSVReader;
import au.com.bytecode.opencsv.CSVWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.regex.Pattern;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import umcg.genetica.collections.intervaltree.NamedGenomicRange;
import umcg.genetica.containers.Pair;

/**
 *
 * @author patri
 */
public class MetaTransCorrelationQtl {

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

		OptionBuilder.withArgName("paths");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Per cohort restult files");
		OptionBuilder.withLongOpt("cohortRes");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("c"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Output file");
		OptionBuilder.withLongOpt("output");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("o"));

	}

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException, IOException {

		System.out.println(HEADER);
		System.out.println();
		System.out.flush(); //flush to make sure header is before errors
		try {
			Thread.sleep(25); //Allows flush to complete
		} catch (InterruptedException ex) {
		}

		final String[] perCohortRestultPaths;
		final String outputPath;

		try {
			final CommandLine commandLine = new PosixParser().parse(OPTIONS, args, false);

			perCohortRestultPaths = commandLine.getOptionValues("c");
			outputPath = commandLine.getOptionValue("o");

		} catch (ParseException ex) {
			System.err.println("Invalid command line arguments: ");
			System.err.println(ex.getMessage());
			System.err.println();
			new HelpFormatter().printHelp(" ", OPTIONS);
			System.exit(1);
			return;
		}

		HashMap<String, MetaTRansCorrelationQtlResult> metaResults = new LinkedHashMap<>();

		for (String cohortResultPath : perCohortRestultPaths) {

			File cohortFile = new File(cohortResultPath);

			CSVReader geneMapReader = new CSVReader(new InputStreamReader(new FileInputStream(cohortFile), ENCODING), '\t', '\0', 1);
			String[] nextLine;
			HashMap<String, ArrayList<NamedGenomicRange>> genes = new HashMap<>();
			while ((nextLine = geneMapReader.readNext()) != null) {
				String snp = nextLine[0];
				double zscoreQtl = Double.parseDouble(nextLine[1]);
				String localGene = nextLine[2];
				String transGene = nextLine[3];

				double zscore = Double.parseDouble(nextLine[5]);
				int samples = Integer.parseInt(nextLine[6]);
				String allele = nextLine[7];

				String key = snp + "-" + localGene + "-" + transGene;

				MetaTRansCorrelationQtlResult result = metaResults.getOrDefault(key, new MetaTRansCorrelationQtlResult(snp, zscoreQtl, localGene, transGene, allele));
				result.addCohortRestult(allele, samples, zscore, cohortFile.getName());
				metaResults.putIfAbsent(key, result);

			}

		}

		CSVWriter outputWriter = new CSVWriter(new FileWriter(outputPath), '\t', CSVWriter.NO_QUOTE_CHARACTER);

		final String[] outputLine = new String[10];
		int c = 0;
		outputLine[c++] = "TransSnp";
		outputLine[c++] = "Z-score";
		outputLine[c++] = "NearbyGene";
		outputLine[c++] = "TargetGene";
		outputLine[c++] = "Cohorts";
		outputLine[c++] = "Samples";
		outputLine[c++] = "Zscores";
		outputLine[c++] = "TotalSamples";
		outputLine[c++] = "MetaZscore";
		outputLine[c++] = "MetaPvalue";
		
		outputWriter.writeNext(outputLine);

		for (MetaTRansCorrelationQtlResult result : metaResults.values()) {

			c = 0;
			outputLine[c++] = result.getTransSnp();
			outputLine[c++] = String.valueOf(result.getZscore());
			outputLine[c++] = result.getLocalGene();
			outputLine[c++] = result.getTransGene();
			outputLine[c++] = result.getCohorts();
			outputLine[c++] = result.getCohortSizes();
			outputLine[c++] = result.getCohortZscores();
			outputLine[c++] = String.valueOf(result.getTotalSampleCount());
			outputLine[c++] = String.valueOf(result.calculateMetaZscore());
			outputLine[c++] = String.valueOf(result.calculateMetaP());

			outputWriter.writeNext(outputLine);

		}
		
		outputWriter.close();

	}

}
