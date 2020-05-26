/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlinteractionanalyser.eqtlinteractionanalyser;

import java.io.*;

import org.apache.commons.cli.*;
import umcg.genetica.io.text.TextFile;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashMap;

/**
 *
 * @author lude
 */
public class EQTLInteractionAnalyser {

    private static final DateFormat DATE_TIME_FORMAT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
    private static final Date currentDataTime = new Date();
    protected static final Options OPTIONS;

    static {

        OPTIONS = new Options();

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("Path to the folder containing expression and genotype data");
        OptionBuilder.withLongOpt("input");
        OptionBuilder.isRequired();
        OPTIONS.addOption(OptionBuilder.create("i"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Path to the output folder");
        OptionBuilder.withLongOpt("output");
        OptionBuilder.isRequired();
        OPTIONS.addOption(OptionBuilder.create("o"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Path to the eQTL file to test for interactions");
        OptionBuilder.withLongOpt("eqtls");
        OPTIONS.addOption(OptionBuilder.create("e"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Path to the eQTL file to correct covariates");
        OptionBuilder.withLongOpt("eqtlsCovariates");
        OPTIONS.addOption(OptionBuilder.create("ec"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Path to the gene annotation file in the format of eQTL mapping pipeline");
        OptionBuilder.withLongOpt("annot");
        OPTIONS.addOption(OptionBuilder.create("a"));

        OptionBuilder.withArgName("int");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Maximum number of covariates to regress out");
        OptionBuilder.withLongOpt("maxcov");
        OPTIONS.addOption(OptionBuilder.create("n"));

        OptionBuilder.withDescription("Interpret the z-score matrices");
        OptionBuilder.withLongOpt("interpret");
        OPTIONS.addOption(OptionBuilder.create("it"));

        OptionBuilder.withDescription("Run permutation");
        OptionBuilder.withLongOpt("permute");
        OPTIONS.addOption(OptionBuilder.create("perm"));

        OptionBuilder.withDescription("Find chi2sum differences for each covariate between 2 consequtive interaction runs");
        OptionBuilder.withLongOpt("chi2sumDiff");
        OPTIONS.addOption(OptionBuilder.create("dif"));

        OptionBuilder.withArgName("int");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Start round for chi2sumDiff option");
        OptionBuilder.withLongOpt("start");
        OPTIONS.addOption(OptionBuilder.create("s"));

        OptionBuilder.withDescription("Preprocess the data");
        OptionBuilder.withLongOpt("preprocess");
        OPTIONS.addOption(OptionBuilder.create("p"));

        OptionBuilder.withDescription("Convert matrix");
        OptionBuilder.withLongOpt("convertMatrix");
        OPTIONS.addOption(OptionBuilder.create("cm"));

        OptionBuilder.withDescription("Skip all normalization step. n must be 1");
        OptionBuilder.withLongOpt("noNormalization");
        OPTIONS.addOption(OptionBuilder.create("nn"));

        OptionBuilder.withDescription("Skip covariate normalization step. n must be 1");
        OptionBuilder.withLongOpt("noCovNormalization");
        OPTIONS.addOption(OptionBuilder.create("ncn"));

        OptionBuilder.withArgName("strings");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("covariates to correct for using an interaction term before running the interaction analysis");
        OptionBuilder.withLongOpt("cov");
        OPTIONS.addOption(OptionBuilder.create("c"));

        OptionBuilder.withArgName("strings");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("Covariates to correct for without interaction term before running the interaction analysis");
        OptionBuilder.withLongOpt("cov2");
        OPTIONS.addOption(OptionBuilder.create("c2"));

        OptionBuilder.withArgName("strings");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("Covariates to correct for without interaction term before running the interaction analysis");
        OptionBuilder.withLongOpt("cohorts");
        OPTIONS.addOption(OptionBuilder.create("ch"));

        OptionBuilder.withArgName("strings");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("Covariates to to test in interaction analysis. Optional, all are tested if not used");
        OptionBuilder.withLongOpt("covTest");
        OPTIONS.addOption(OptionBuilder.create("ct"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("File containing the covariates to correct for using an interaction term before running the interaction analysis. No header, each covariate on a separate line");
        OptionBuilder.withLongOpt("covFile");
        OPTIONS.addOption(OptionBuilder.create("cf"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("File containing the SNPs to swap");
        OptionBuilder.withLongOpt("swap");
        OPTIONS.addOption(OptionBuilder.create("sw"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Included samples");
        OptionBuilder.withLongOpt("includedSamples");
        OPTIONS.addOption(OptionBuilder.create("is"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Gene annotation file");
        OptionBuilder.withLongOpt("geneAnnotation");
        OPTIONS.addOption(OptionBuilder.create("ga"));

        OptionBuilder.withArgName("int");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Number of threads");
        OptionBuilder.withLongOpt("threads");
        OPTIONS.addOption(OptionBuilder.create("nt"));

        OptionBuilder.withArgName("int");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Z-score difference threshold for interpretation");
        OptionBuilder.withLongOpt("threshold");
        OPTIONS.addOption(OptionBuilder.create("thr"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("SNPs to test");
        OptionBuilder.withLongOpt("snpsToTest");
        OPTIONS.addOption(OptionBuilder.create("snps"));

        OptionBuilder.withArgName("int");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Number of PCs to correct for");
        OptionBuilder.withLongOpt("numpc");
        OPTIONS.addOption(OptionBuilder.create("pc"));
    }

    public static void main(String[] args) throws IOException, Exception {
        System.out.println("Starting interaction analysis");
        System.out.println("Current date and time: " + DATE_TIME_FORMAT.format(currentDataTime));
        System.out.println();

        String inputDir, outputDir, eqtlFile = null, annotationFile = null;
        final File snpsToSwapFile;
        int maxNumCovariatesToRegress = 20;
        int numPCsToRegress = 25;
        int numThreads;
        final boolean interpret, chi2sumDiff, permute, preproces;
        final int startRoundCompareChi2, threshold;

        HashMap hashSamples;

        final String[] covariates;
        final String[] covariates2;
        final String[] cohorts;
        final String[] covariatesToTest;
        final File ensgAnnotationFile;
        final File snpsToTestFile;
        final boolean skipNormalization;
        final boolean skipCovariateNormalization;
        final boolean convertMatrix;
        final String eqtlFileCovariates;

        try {
            final CommandLine commandLine = new PosixParser().parse(OPTIONS, args, false);

            inputDir = commandLine.getOptionValue("i");
            outputDir = commandLine.getOptionValue("o");

            if (commandLine.hasOption('e')) {
                eqtlFile = commandLine.getOptionValue("e");
            }

            eqtlFileCovariates = commandLine.getOptionValue("ec", null);

            if (commandLine.hasOption('n')) {
                maxNumCovariatesToRegress = Integer.parseInt(commandLine.getOptionValue("n"));
            }
            if (commandLine.hasOption("thr")) {
                threshold = Integer.parseInt(commandLine.getOptionValue("thr"));
            } else {
                threshold = 3;
            }

            interpret = commandLine.hasOption("it");
            chi2sumDiff = commandLine.hasOption("dif");
            permute = commandLine.hasOption("perm");
            preproces = commandLine.hasOption("p");
            convertMatrix = commandLine.hasOption("cm");

            if (commandLine.hasOption('s')) {
                startRoundCompareChi2 = Integer.parseInt(commandLine.getOptionValue("s"));
            } else if (chi2sumDiff) {
                throw new Exception("Set -s");
            } else {
                startRoundCompareChi2 = 0;
            }

            if (commandLine.hasOption('a')) {
                annotationFile = commandLine.getOptionValue("a");
            }

            if (commandLine.hasOption("cf")) {
                TextFile covFile = new TextFile(commandLine.getOptionValue("cf"), false);
                covariates = covFile.readAsArray();
                covFile.close();
            } else if (commandLine.hasOption("c")) {
                covariates = commandLine.getOptionValues("c");
            } else {
                covariates = new String[0];
            }

            if (commandLine.hasOption("c2")) {
                covariates2 = commandLine.getOptionValues("c2");
            } else {
                covariates2 = new String[0];
            }

            if (commandLine.hasOption("ch")) {
                cohorts = commandLine.getOptionValues("ch");
            } else {
                cohorts = null;
            }

            if (commandLine.hasOption("ct")) {
                covariatesToTest = commandLine.getOptionValues("ct");
            } else {
                covariatesToTest = null;
            }

            if (commandLine.hasOption("sw")) {
                snpsToSwapFile = new File(commandLine.getOptionValue("sw"));
            } else {
                snpsToSwapFile = null;
            }

            if (commandLine.hasOption("snps")) {
                snpsToTestFile = new File(commandLine.getOptionValue("snps"));
            } else {
                snpsToTestFile = null;
            }

            skipNormalization = commandLine.hasOption("nn");
            if (skipNormalization && maxNumCovariatesToRegress != 1) {
                System.err.println("n must be one if normalization is turned off");
                System.exit(-1);
            }

            skipCovariateNormalization = commandLine.hasOption("ncn");
            if (skipCovariateNormalization && maxNumCovariatesToRegress != 1) {
                System.err.println("n must be one if covariate normalization is turned off");
                System.exit(-1);
            }

            if (commandLine.hasOption("is")) {
                File samplesToIncludeFile = new File(commandLine.getOptionValue("is"));
                System.out.println("Samples to include file: " + samplesToIncludeFile.getAbsolutePath());
                hashSamples = new HashMap();
                BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(samplesToIncludeFile), "UTF-8"));
                String line;
                while ((line = reader.readLine()) != null) {
                    hashSamples.put(line, null);
                    hashSamples.put(line + "_exp", null);
                    hashSamples.put(line + "_dosage", null);
                }
            } else {
                hashSamples = null;
            }

            if (commandLine.hasOption("ga")) {
                ensgAnnotationFile = new File(commandLine.getOptionValue("ga"));
            } else {
                ensgAnnotationFile = null;
            }
            if (commandLine.hasOption("nt")) {
                numThreads = Integer.parseInt(commandLine.getOptionValue("nt"));
            } else {
                numThreads = Runtime.getRuntime().availableProcessors();
            }
            if (commandLine.hasOption("pc")) {
                numPCsToRegress = Integer.parseInt(commandLine.getOptionValue("pc"));
            }

        } catch (ParseException ex) {
            System.err.println("Invalid command line arguments: ");
            System.err.println(ex.getMessage());
            System.err.println();
            new HelpFormatter().printHelp(" ", OPTIONS);
            System.exit(1);
            return;
        }

        if (preproces) {
            TestEQTLDatasetForInteractions interactor = new TestEQTLDatasetForInteractions(inputDir, outputDir);
            interactor.preprocessData();
        } else if (interpret) {
            TestEQTLDatasetForInteractions interactor = new TestEQTLDatasetForInteractions(inputDir, outputDir);
            interactor.interpretInteractionZScoreMatrix(maxNumCovariatesToRegress, startRoundCompareChi2, threshold);
        } else if (chi2sumDiff) {
            TestEQTLDatasetForInteractions interactor = new TestEQTLDatasetForInteractions(inputDir, outputDir);
            interactor.findChi2SumDifferences(maxNumCovariatesToRegress, startRoundCompareChi2, ensgAnnotationFile);
        } else if (convertMatrix) {
            System.out.println("input file: " + inputDir);
            System.out.println("output file: " + outputDir);
            if (inputDir.equals(outputDir)) {
                System.err.println("Error, input == output");
                System.exit(1);
            }
            new ExpressionDataset(inputDir).save(outputDir);
        } else {
            new TestEQTLDatasetForInteractions(inputDir, outputDir, eqtlFile, maxNumCovariatesToRegress, annotationFile, covariates, covariates2, snpsToSwapFile, permute, covariatesToTest, hashSamples, numThreads, cohorts, snpsToTestFile, skipNormalization, skipCovariateNormalization, eqtlFileCovariates, numPCsToRegress);
        }
    }

}
