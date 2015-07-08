/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package nl.systemsgenetics.eqtlinteractionanalyser.eqtlinteractionanalyser;

import com.sun.org.apache.xpath.internal.operations.Bool;
import org.apache.commons.cli.*;

import java.io.File;
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

        Option option;

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

        OptionBuilder.withArgName("int");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Maximum number of covariates to regress out");
        OptionBuilder.withLongOpt("maxcov");
        OPTIONS.addOption(OptionBuilder.create("n"));

        OptionBuilder.withArgName("boolean");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Interpret the z-score matrices");
        OptionBuilder.withLongOpt("interpret");
        OPTIONS.addOption(OptionBuilder.create("t"));
    }

    public static void main(String[] args) throws IOException {
        System.out.println("Starting interaction analysis");
        System.out.println("Current date and time: " + DATE_TIME_FORMAT.format(currentDataTime));
        System.out.println();


        System.out.flush(); //flush to make sure header is before errors
        try {
            Thread.sleep(25); //Allows flush to complete
        } catch (InterruptedException ex) {
        }

        String inputDir, outputDir, eqtlFile = null;
        int maxNumCovariatesToRegress = 300;
        boolean interpret = false;
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
            if (commandLine.hasOption('t')) {
                interpret = Boolean.parseBoolean(commandLine.getOptionValue("t"));
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
            interactor.interpretInteractionZScoreMatrix();
        }
        else {
            new TestEQTLDatasetForInteractions(inputDir, outputDir, eqtlFile, maxNumCovariatesToRegress);
        }
    }
    
}
