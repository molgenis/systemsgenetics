/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package nl.systemsgenetics.eqtlinteractionanalyser.eqtlinteractionanalyser;

import org.apache.commons.cli.*;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;

/**
 *
 * @author lude
 */
public class EQTLInteractionAnalyser {

    private static final DateFormat DATE_TIME_FORMAT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
    private static final Date currentDataTime = new Date();
    private static final Options OPTIONS;

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
        OptionBuilder.withDescription("Path to the gene annotation file in the format of eQTL mapping pipeline");
        OptionBuilder.withLongOpt("annot");
        OPTIONS.addOption(OptionBuilder.create("a"));

        OptionBuilder.withArgName("int");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Maximum number of covariates to regress out");
        OptionBuilder.withLongOpt("maxcov");
        OPTIONS.addOption(OptionBuilder.create("n"));

        OptionBuilder.withArgName("boolean");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Interpret the z-score matrices");
        OptionBuilder.withLongOpt("interpret");
        OPTIONS.addOption(OptionBuilder.create("it"));

        OptionBuilder.withArgName("boolean");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Find chi2sum differences for each covariate between 2 consequtive interaction runs");
        OptionBuilder.withLongOpt("chi2sumDiff");
        OPTIONS.addOption(OptionBuilder.create("dif"));

        OptionBuilder.withArgName("string");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("covariates to correct for before running the interaction analysis");
        OptionBuilder.withLongOpt("cov");
        OPTIONS.addOption(OptionBuilder.create("c"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("File containing the covariates to correct for before running the interaction analysis. No header, each covariate on a separate line");
        OptionBuilder.withLongOpt("covFile");
        OPTIONS.addOption(OptionBuilder.create("cf"));
    }

    public static void main(String[] args) throws IOException {
        System.out.println("Starting interaction analysis");
        System.out.println("Current date and time: " + DATE_TIME_FORMAT.format(currentDataTime));
        System.out.println();
        
        String inputDir, outputDir, eqtlFile = null, annotationFile = null;
        int maxNumCovariatesToRegress = 20;
        boolean interpret = false, chi2sumDiff = false;
        String[] covariates = null;
        try {
            final CommandLine commandLine = new PosixParser().parse(OPTIONS, args, false);

            inputDir = commandLine.getOptionValue("i");
            outputDir = commandLine.getOptionValue("o");

            if (commandLine.hasOption('e')) {
                eqtlFile = commandLine.getOptionValue("e");
            }
            if (commandLine.hasOption('n')) {
                maxNumCovariatesToRegress = Integer.parseInt(commandLine.getOptionValue("n"));
            }
            if (commandLine.hasOption("it")) {
                interpret = Boolean.parseBoolean(commandLine.getOptionValue("t"));
            }
            if (commandLine.hasOption("dif")) {
                chi2sumDiff = Boolean.parseBoolean(commandLine.getOptionValue("d"));
            }
            if (commandLine.hasOption('a')) {
                annotationFile = commandLine.getOptionValue("a");
            }

            if (commandLine.hasOption("cf")) {
                TextFile covFile = new TextFile(commandLine.getOptionValue("cf"), false);
                covariates = covFile.readAsArray();
                covFile.close();
            }
            else if (commandLine.hasOption("c")){
                covariates = commandLine.getOptionValues("cf");
            }

        } catch (ParseException ex) {
            System.err.println("Invalid command line arguments: ");
            System.err.println(ex.getMessage());
            System.err.println();
            new HelpFormatter().printHelp(" ", OPTIONS);
            System.exit(1);
            return;
        }

        if (interpret){
            TestEQTLDatasetForInteractions interactor = new TestEQTLDatasetForInteractions(inputDir, outputDir);
            interactor.interpretInteractionZScoreMatrix(maxNumCovariatesToRegress);
        }
        else if (chi2sumDiff){
            TestEQTLDatasetForInteractions interactor = new TestEQTLDatasetForInteractions(inputDir, outputDir);
            interactor.findChi2SumDifferences(maxNumCovariatesToRegress);
        }
        else {
            new TestEQTLDatasetForInteractions(inputDir, outputDir, eqtlFile, maxNumCovariatesToRegress, annotationFile, covariates);
        }
    }

}
