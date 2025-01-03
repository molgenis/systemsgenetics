/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.runners.options;

import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.zip.GZIPInputStream;
import nl.systemsgenetics.downstreamer.gene.Gene;
import nl.systemsgenetics.downstreamer.pathway.PathwayDatabase;
import static nl.systemsgenetics.downstreamer.runners.options.OptionsBase.OPTIONS;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.mahout.math.Arrays;

/**
 *
 * @author patri
 */
public class OptionsModeEnrichment extends OptionsBase {

	private static final Logger LOGGER = LogManager.getLogger(OptionsBase.class);

	private final File covariates;
	private final List<PathwayDatabase> pathwayDatabases;
	private final boolean forceNormalGenePvalues;
	private final boolean forceNormalPathwayPvalues;
	private final boolean regressGeneLengths;
	private final File geneInfoFile;
	private final File singleGwasFile;
	private final String gwasPvalueMatrixPath;
	private final boolean excludeHla;
	private final boolean skipPvalueToZscore;
	private final String geneGeneCorrelationPrefix;
	private final boolean unitTestMode;//only set true for unit testing
	private final double fdrThresholdEigenvectors;

	static {

		OptionBuilder.withArgName("boolean");
		OptionBuilder.withDescription("Linearly correct for the effect of gene lengths on gene p-values");
		OptionBuilder.withLongOpt("regress-gene-lengths");
		OPTIONS.addOption(OptionBuilder.create("rgl"));

		OptionBuilder.withArgName("boolean");
		OptionBuilder.withDescription("Force normal gene p-values before pathway enrichtment");
		OptionBuilder.withLongOpt("forceNormalGenePvalues");
		OPTIONS.addOption(OptionBuilder.create("fngp"));

		OptionBuilder.withArgName("boolean");
		OptionBuilder.withDescription("Force normal pathway scores / eigen vectors before pathway enrichtment");
		OptionBuilder.withLongOpt("forceNormalPathwayPvalues");
		OPTIONS.addOption(OptionBuilder.create("fnpp"));

		OptionBuilder.withArgName("name=path");
		OptionBuilder.hasArgs();
		OptionBuilder.withValueSeparator();
		OptionBuilder.withDescription("Pathway databases, .dat or .datg matrix with either z-scores for predicted gene pathway associations or 0 / 1 for gene assignments");
		OptionBuilder.withLongOpt("pathwayDatabase");
		OPTIONS.addOption(OptionBuilder.create("pd"));

		OptionBuilder.withArgName("name=path");
		OptionBuilder.hasArgs();
		OptionBuilder.withValueSeparator();
		OptionBuilder.withDescription("Expression eigen vectors for gene prioritizaion. In .dat or .datg foramat");
		OptionBuilder.withLongOpt("expressionEigenVectors");
		OPTIONS.addOption(OptionBuilder.create("eev"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Instead of name=path for pathways or expession eigenvectors use file with 3 columns (no header): name<tab>path<tab>true/false. Use true in the last column to indicate if eigenvectors instead of pathways.");
		OptionBuilder.withLongOpt("pathwayEigenFile");
		OPTIONS.addOption(OptionBuilder.create("pef"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("File with gene info. col1: geneName (ensg) col2: chr col3: startPos col4: stopPos col5: geneType col6: chrArm");
		OptionBuilder.withLongOpt("genes");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("ge"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Path to files with gene-gene corelations. Specify character until 'chr_arm'.");
		OptionBuilder.withLongOpt("geneCorrelations");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("gc"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("File with covariates used to correct the gene p-values. Works in conjunction with -rgl. Residuals of this regression are used as input for the GLS");
		OptionBuilder.withLongOpt("covariates");
		OPTIONS.addOption(OptionBuilder.create("cov"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("GWAS gene p-values. Rows genes, Cols: genes, pvalue, nSNPs, min SNP p-value. Tab seperated txt or txt.gz file or .dat / datg binary matrix format");
		OptionBuilder.withLongOpt("gwas");
		OPTIONS.addOption(OptionBuilder.create("g"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Expects 3 files with equal and same ordered columns (traits) and rows (genes). <path>_pvalues with gene p-values; <path>_nvar number of variants per gene; <path>_minVarPvalue. Each file is Tab seperated txt/txt.gz file or .dat / datg binary matrix format. The number of variants per gene is assumed to be equal, it is allowed ot have only a single column in this file.");
		OptionBuilder.withLongOpt("gwasMatrix");
		OPTIONS.addOption(OptionBuilder.create("gm"));

		OptionBuilder.withArgName("boolean");
		OptionBuilder.withDescription("Exclude HLA locus during pathway enrichments (chr6 20mb - 40mb)");
		OptionBuilder.withLongOpt("excludeHla");
		OPTIONS.addOption(OptionBuilder.create("eh"));

		OptionBuilder.withArgName("boolean");
		OptionBuilder.withDescription("Set true to skip converting gene p-values to z-scores. Only do this if the gene scores are already z-scores");
		OptionBuilder.withLongOpt("noPvalueToZscore");
		OPTIONS.addOption(OptionBuilder.create("nptz"));

		OptionBuilder.withArgName("double");
		OptionBuilder.withDescription("FDR threshold to determine which eigenvectors should be used for key-gene scores");
		OptionBuilder.withLongOpt("eigenFdr");
		OptionBuilder.hasArg();
		OPTIONS.addOption(OptionBuilder.create("efdr"));

	}

	/**
	 * Constructor used for the unit testing
	 *
	 * @param covariates
	 * @param pathwayDatabases
	 * @param forceNormalGenePvalues
	 * @param forceNormalPathwayPvalues
	 * @param regressGeneLengths
	 * @param geneInfoFile
	 * @param singleGwasFile
	 * @param gwasPvalueMatrixPath
	 * @param excludeHla
	 * @param skipPvalueToZscore
	 * @param geneGeneCorrelationPrefix
	 * @param numberOfThreadsToUse
	 * @param outputBasePath
	 * @param logFile
	 * @param mode
	 * @param debugMode
	 * @param jblas
	 * @param unitTestMode
	 */
	public OptionsModeEnrichment(File covariates, List<PathwayDatabase> pathwayDatabases, boolean forceNormalGenePvalues, boolean forceNormalPathwayPvalues, boolean regressGeneLengths, File geneInfoFile, File singleGwasFile, String gwasPvalueMatrixPath, boolean excludeHla, boolean skipPvalueToZscore, String geneGeneCorrelationPrefix, int numberOfThreadsToUse, File outputBasePath, File logFile, DownstreamerMode mode, boolean debugMode, boolean jblas, boolean unitTestMode, double fdrThresholdEigenvectors) {
		super(numberOfThreadsToUse, outputBasePath, logFile, mode, debugMode, jblas);
		this.covariates = covariates;
		this.pathwayDatabases = pathwayDatabases;
		this.forceNormalGenePvalues = forceNormalGenePvalues;
		this.forceNormalPathwayPvalues = forceNormalPathwayPvalues;
		this.regressGeneLengths = regressGeneLengths;
		this.geneInfoFile = geneInfoFile;
		this.singleGwasFile = singleGwasFile;
		this.gwasPvalueMatrixPath = gwasPvalueMatrixPath;
		this.excludeHla = excludeHla;
		this.skipPvalueToZscore = skipPvalueToZscore;
		this.geneGeneCorrelationPrefix = geneGeneCorrelationPrefix;
		this.unitTestMode = unitTestMode;
		this.fdrThresholdEigenvectors = fdrThresholdEigenvectors;
	}

	public OptionsModeEnrichment(String[] args) throws ParseException, IOException {
		super(args);

		unitTestMode = false;

		// Parse arguments
		final CommandLineParser parser = new PosixParser();
		final CommandLine commandLine = parser.parse(OPTIONS, args, false);

		if (commandLine.hasOption("cov")) {
			covariates = new File(commandLine.getOptionValue("cov"));
		} else {
			covariates = null;
		}

		pathwayDatabases = parsePathwaysAndExpressionEigenVectors(commandLine);
		forceNormalGenePvalues = commandLine.hasOption("fngp");
		forceNormalPathwayPvalues = commandLine.hasOption("fnpp");
		regressGeneLengths = commandLine.hasOption("rgl");
		skipPvalueToZscore = commandLine.hasOption("nptz");
		geneInfoFile = new File(commandLine.getOptionValue("ge"));
		geneGeneCorrelationPrefix = commandLine.getOptionValue("gc");

		if (commandLine.hasOption("g") && commandLine.hasOption("gm")) {
			throw new ParseException("Provide either -g or -gm but not both");
		} else if (commandLine.hasOption("g")) {
			singleGwasFile = new File(commandLine.getOptionValue('g'));
			gwasPvalueMatrixPath = null;
		} else if (commandLine.hasOption("gm")) {
			singleGwasFile = null;
			gwasPvalueMatrixPath = commandLine.getOptionValue("gm");
		} else {
			throw new ParseException("Provide either -g or -gm");
		}

		excludeHla = commandLine.hasOption("eh");
		
		if(commandLine.hasOption("eigenFdr")){
			fdrThresholdEigenvectors = Double.parseDouble(commandLine.getOptionValue("eigenFdr"));
		} else {
			fdrThresholdEigenvectors = 0.05;
		}

		printOptions();

	}

	private static List<PathwayDatabase> parsePathwaysAndExpressionEigenVectors(final CommandLine commandLine) throws ParseException, FileNotFoundException, IOException {

		final List<PathwayDatabase> pathwayDatabases = new ArrayList<>();
		final HashSet<String> duplicateChecker = new HashSet<>();

		if (commandLine.hasOption("pd")) {

			String[] pdValues = commandLine.getOptionValues("pd");

			if (pdValues.length % 2 != 0) {
				throw new ParseException("Error parsing --pathwayDatabase. Must be in name=path format");
			}

			for (int i = 0; i < pdValues.length; i += 2) {

				if (!duplicateChecker.add(pdValues[i])) {
					throw new ParseException("Error parsing --pathwayDatabase. Duplicate database name found");
				}

				pathwayDatabases.add(new PathwayDatabase(pdValues[i], pdValues[i + 1], false));

			}

		}

		if (commandLine.hasOption("eev")) {

			String[] pdValues = commandLine.getOptionValues("eev");

			if (pdValues.length % 2 != 0) {
				throw new ParseException("Error parsing --ExpressionEigenVectors. Must be in name=path format");
			}

			for (int i = 0; i < pdValues.length; i += 2) {

				if (!duplicateChecker.add(pdValues[i])) {
					throw new ParseException("Error parsing --ExpressionEigenVectors. Duplicate database name found, note must also be unique with --pathwayDatabase");
				}

				pathwayDatabases.add(new PathwayDatabase(pdValues[i], pdValues[i + 1], true));

			}

		}

		if (commandLine.hasOption("pef")) {

			File pathwayFile = new File(commandLine.getOptionValue("pef"));

			final CSVReader reader
					= new CSVReaderBuilder((new BufferedReader(
							pathwayFile.getName().endsWith(".gz")
							? new InputStreamReader(new GZIPInputStream(new FileInputStream(pathwayFile)))
							: new FileReader(pathwayFile)
					))).withSkipLines(0).withCSVParser(
							new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build()
					).build();

			String[] nextLine;
			while ((nextLine = reader.readNext()) != null) {

				if (!duplicateChecker.add(nextLine[0])) {
					throw new ParseException("Error parsing --pathwayEigenFile. Duplicate database name found, note must also be unique with --pathwayDatabase --ExpressionEigenVectors ");
				}
				
				if(nextLine.length != 3){
					throw new ParseException("Error parsing --pathwayEigenFile. Each row should contain 3 columns but found: " + nextLine.length + " at line: " + Arrays.toString(nextLine));
				}

				pathwayDatabases.add(new PathwayDatabase(nextLine[0], nextLine[1], Boolean.parseBoolean(nextLine[2])));

			}

		}

		if (duplicateChecker.isEmpty()) {

			throw new ParseException("Mode ENRICH requirers either --ExpressionEigenVectors or --pathwayDatabase");

		}

		return Collections.unmodifiableList(pathwayDatabases);
	}

	@Override
	public void printOptions() {
		super.printOptions();

		if (singleGwasFile != null) {
			LOGGER.info(" * GWAS gene p-value, variant count, min variant p-value file: " + singleGwasFile.getPath());
		} else {
			LOGGER.info(" * GWAS gene p-value matrix path: " + gwasPvalueMatrixPath);
		}

		if (covariates != null) {
			LOGGER.info(" * Covariates: " + covariates.getPath());
		}

		LOGGER.info(" * Pathway databases: ");
		for (PathwayDatabase curDb : pathwayDatabases) {
			if (!curDb.isEigenvectors()) {
				LOGGER.info("    - " + curDb.getName() + "\t" + curDb.getLocation());
			}
		}
		LOGGER.info(" * Expression eigenvectors: ");
		for (PathwayDatabase curDb : pathwayDatabases) {
			if (curDb.isEigenvectors()) {
				LOGGER.info("    - " + curDb.getName() + "\t" + curDb.getLocation());
			}
		}

		LOGGER.info(" * Gene info file: " + geneInfoFile.getPath());
		LOGGER.info(" * Path to gene-gene correlation files: " + geneGeneCorrelationPrefix);
		LOGGER.info(" * Do inverse force normal of gene p-values: " + forceNormalGenePvalues);
		LOGGER.info(" * Do inverse force normal of pathway scores: " + forceNormalPathwayPvalues);
		LOGGER.info(" * Skip gene p-values to Z-score (should only be disabled if already z-score per gene): " + skipPvalueToZscore);
		LOGGER.info(" * Correct gene p-values for gene length: " + regressGeneLengths);
		LOGGER.info(" * Exclude HLA during enrichment analysis: " + (excludeHla ? "on" : "off"));
		LOGGER.info(" * FDR threshold for eigenvectors to use: " + fdrThresholdEigenvectors);

	}

	public boolean isForceNormalGenePvalues() {
		return forceNormalGenePvalues;
	}

	public List<PathwayDatabase> getPathwayDatabases() {
		return pathwayDatabases;
	}

	public File getCovariates() {
		return covariates;
	}

	public boolean isForceNormalPathwayPvalues() {
		return forceNormalPathwayPvalues;
	}

	public boolean isRegressGeneLengths() {
		return regressGeneLengths;
	}

	public File getGeneInfoFile() {
		return geneInfoFile;
	}

	public File getSingleGwasFile() {
		return singleGwasFile;
	}

	public String getGwasPvalueMatrixPath() {
		return gwasPvalueMatrixPath;
	}

	public boolean isExcludeHla() {
		return excludeHla;
	}

	public boolean isSkipPvalueToZscore() {
		return skipPvalueToZscore;
	}

	public String getGeneGeneCorrelationPrefix() {
		return geneGeneCorrelationPrefix;
	}

	public boolean isUnitTestMode() {
		return unitTestMode;
	}

	public double getFdrThresholdEigenvectors() {
		return fdrThresholdEigenvectors;
	}

}
