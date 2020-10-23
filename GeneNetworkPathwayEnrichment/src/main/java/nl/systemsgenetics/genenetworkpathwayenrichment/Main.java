package nl.systemsgenetics.genenetworkpathwayenrichment;

import org.apache.commons.cli.*;
import umcg.genetica.io.Gpio;

import java.io.IOException;

public class Main {

    private static final Options OPTIONS;

    public static void main(String[] args) {
        CommandLineParser parser = new PosixParser();
        try {
            final CommandLine commandLine = parser.parse(OPTIONS, args, false);

            GeneSetEnrichmentWilcoxon gsea = new GeneSetEnrichmentWilcoxon();

            String matrix = null;
            String genesetfile = null;
            String annotationfile = null;
            String ensgtohugofile = null;
            String limitgenes = null;
            String meanexpfile = null;
            boolean randomlymatchnumberofgenes = false;
            String output = null;

            boolean run = true;
            matrix = commandLine.getOptionValue("m");
            if (Gpio.exists(matrix)) {
                System.out.println(matrix + " does not exist.");
                run = false;
            }
            genesetfile = commandLine.getOptionValue("g");
            if (Gpio.exists(genesetfile)) {
                System.out.println(genesetfile + " does not exist.");
                run = false;
            }
            if (commandLine.hasOption("p")) {
                annotationfile = commandLine.getOptionValue("p");
                if (!Gpio.exists(annotationfile)) {
                    System.out.println(annotationfile + " does not exist.");
                    run = false;
                }
            }
            if (commandLine.hasOption("a")) {
                ensgtohugofile = commandLine.getOptionValue("a");
                if (!Gpio.exists(annotationfile)) {
                    System.out.println(ensgtohugofile + " does not exist.");
                    run = false;
                }
            }
            if (commandLine.hasOption("l")) {
                limitgenes = commandLine.getOptionValue("l");
                if (!Gpio.exists(limitgenes)) {
                    System.out.println(limitgenes + " does not exist.");
                    run = false;
                }
            }
            if (commandLine.hasOption("e")) {
                meanexpfile = commandLine.getOptionValue("e");
                if (!Gpio.exists(meanexpfile)) {
                    System.out.println(meanexpfile + " does not exist.");
                    run = false;
                }
            }

            randomlymatchnumberofgenes = commandLine.hasOption("r");
            output = commandLine.getOptionValue("o");
            if (run) {
                gsea.run(matrix, genesetfile, annotationfile, ensgtohugofile, limitgenes, meanexpfile, randomlymatchnumberofgenes, output);
            } else {
                printHelp();
            }
        } catch (ParseException e) {
            System.out.println(e.getMessage());
            printHelp();
        } catch (IOException e) {
            e.printStackTrace();
        }


    }


    static {

        OPTIONS = new Options();

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("Pathway z-score matrix");
        OptionBuilder.withLongOpt("matrix");
        OptionBuilder.isRequired();
        OPTIONS.addOption(OptionBuilder.create("m"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("File listing (foreground) genes to test");
        OptionBuilder.withLongOpt("genes");
        OptionBuilder.isRequired();
        OPTIONS.addOption(OptionBuilder.create("g"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("File with annotations per ontology/pathway term");
        OptionBuilder.withLongOpt("pathwayannotation");
        OPTIONS.addOption(OptionBuilder.create("p"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("File with hugo annotations per ENSG gene id");
        OptionBuilder.withLongOpt("geneannotation");
        OPTIONS.addOption(OptionBuilder.create("a"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("File with list of genes on which to limit the analysis");
        OptionBuilder.withLongOpt("genelimit");
        OPTIONS.addOption(OptionBuilder.create("l"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("File with mean expression per gene");
        OptionBuilder.withLongOpt("geneexp");
        OPTIONS.addOption(OptionBuilder.create("e"));

        OptionBuilder.withArgName("path");
        OptionBuilder.withDescription("Randomly select genes as background (does not work with -e option)");
        OptionBuilder.withLongOpt("randombackground");
        OPTIONS.addOption(OptionBuilder.create("r"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("Output path");
        OptionBuilder.withLongOpt("output");
        OptionBuilder.isRequired();
        OPTIONS.addOption(OptionBuilder.create("o"));

    }

    public static void printHelp() {
        System.out.println("Gene set enrichment using GeneNetwork pathway Z-score matrices");
        System.out.println("Different options for background selection exist:\n" +
                "When using -l, foreground and background genes are limited to the set of genes specified. This option also applies to -r and -e.\n" +
                "When using -r, a random set of background genes is selected with equal size to genes specified with -g\n" +
                "When using -e, an equal number of background genes as specified with -g are selected by matching on expression values");
        System.out.println();
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp(" ", OPTIONS);
    }
}
