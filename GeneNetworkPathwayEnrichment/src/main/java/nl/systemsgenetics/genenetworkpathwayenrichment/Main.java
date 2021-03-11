package nl.systemsgenetics.genenetworkpathwayenrichment;

import org.apache.commons.cli.*;
import umcg.genetica.io.Gpio;

import java.io.IOException;

public class Main {

    private static final Options OPTIONS;

    public static void main(String[] args) {
        System.out.println("Gene set enrichment using GeneNetwork pathway Z-score or Identity matrices");
        System.out.println();
        CommandLineParser parser = new PosixParser();
        try {
            final CommandLine commandLine = parser.parse(OPTIONS, args, false);

            String matrix = null;
            String foregroundFile = null;
            String backgroundFile = null;

            String pathwayAnnotationFile = null;
            String ensgtohugofile = null;
            String meanexpfile = null;

            int permutations = 0;
            String output = null;

            boolean run = true;
            matrix = commandLine.getOptionValue("m");
            if (!Gpio.exists(matrix)) {
                System.out.println(matrix + " does not exist.");
                run = false;
            }
            foregroundFile = commandLine.getOptionValue("f");
            if (!Gpio.exists(foregroundFile)) {
                System.out.println(foregroundFile + " does not exist.");
                run = false;
            }

            if (commandLine.hasOption("b")) {
                backgroundFile = commandLine.getOptionValue("b");
                if (!Gpio.exists(backgroundFile)) {
                    System.out.println(backgroundFile + " does not exist.");
                    run = false;
                }
            }

            if (commandLine.hasOption("p")) {
                pathwayAnnotationFile = commandLine.getOptionValue("p");
                if (!Gpio.exists(pathwayAnnotationFile)) {
                    System.out.println(pathwayAnnotationFile + " does not exist.");
                    run = false;
                }
            }

            if (commandLine.hasOption("a")) {
                ensgtohugofile = commandLine.getOptionValue("a");
                if (!Gpio.exists(ensgtohugofile)) {
                    System.out.println(ensgtohugofile + " does not exist.");
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

            if (commandLine.hasOption("perm")) {
                String permutationStr = commandLine.getOptionValue("perm");
                permutations = Integer.parseInt(permutationStr);
            }

            GeneSetEnrichment.TESTTYPE testtype = GeneSetEnrichment.TESTTYPE.WILCOXON;
            if (commandLine.hasOption("type")) {
                String typestr = commandLine.getOptionValue("type");
                if (typestr.equals("wilcox")) {
                    testtype = GeneSetEnrichment.TESTTYPE.WILCOXON;
                } else if (typestr.equals("fet")) {
                    testtype = GeneSetEnrichment.TESTTYPE.FISHEREXACT;
                }
            }


//            if (commandLine.hasOption("x")) {
//                squarevalues = true;
//            }
//
//            if (commandLine.hasOption("lrt")) {
//                lrt = true;
//            }
//
//            if (commandLine.hasOption("properpermute")) {
//                runProperPermutation = true;
//                if (!permute) {
//                    System.out.println("-z/--properpermute requires -p option");
//                    run = false;
//                }
//                if (meanexpfile == null) {
//                    System.out.println("-z/--properpermute requires -e option");
//                    run = false;
//                }
//            }

            boolean randomlymatchnumberofgenes = commandLine.hasOption("r");
            output = commandLine.getOptionValue("o");
            if (run) {

                GeneSetEnrichment g = new GeneSetEnrichment();
                g.setMakeBackgroundEqualSize(randomlymatchnumberofgenes);
                g.run(matrix, foregroundFile, backgroundFile, pathwayAnnotationFile, ensgtohugofile, meanexpfile, permutations, testtype, output);

//                if (lrt) {
//                    GeneSetEnrichmentLogisticRegression lr = new GeneSetEnrichmentLogisticRegression();
//                    lr.run(matrix, genesetfile, limitgenes, annotationfile, ensgtohugofile, meanexpfile, permutations, output);
//                } else
////                if (runProperPermutation) {
//                    gsea.runProperBackgroundCorrection(matrix, genesetfile, annotationfile, ensgtohugofile, limitgenes, meanexpfile, permutations, squarevalues, output);
//
////                } else {
////                    System.out.println("Does not work.");
////                    gsea.run(matrix, genesetfile, annotationfile, ensgtohugofile, limitgenes, meanexpfile, randomlymatchnumberofgenes, permute, permutations, output);
////                }

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
        OptionBuilder.withLongOpt("foreground");
        OptionBuilder.isRequired();
        OPTIONS.addOption(OptionBuilder.create("f"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("File listing background genes");
        OptionBuilder.withLongOpt("background");
        OptionBuilder.isRequired();
        OPTIONS.addOption(OptionBuilder.create("b"));

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
        OptionBuilder.withDescription("File with mean expression per gene");
        OptionBuilder.withLongOpt("meanexp");
        OPTIONS.addOption(OptionBuilder.create("e"));

        OptionBuilder.withArgName("int");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("Number of permutations to run");
        OptionBuilder.withLongOpt("perm");
        OPTIONS.addOption(OptionBuilder.create("perm"));

        OptionBuilder.withArgName("string");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("Test type to run [wilcox|fet]: Wilcoxon or Fisher Exact Test; defaults to Wilcoxon test.");
        OptionBuilder.withLongOpt("type");
        OPTIONS.addOption(OptionBuilder.create("t"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("Output path");
        OptionBuilder.withLongOpt("output");
        OptionBuilder.isRequired();
        OPTIONS.addOption(OptionBuilder.create("o"));

        OptionBuilder.withArgName("path");
        OptionBuilder.withDescription("Randomly select genes as background to match size of foreground.");
        OptionBuilder.withLongOpt("randombackground");
        OPTIONS.addOption(OptionBuilder.create("r"));

        // additional options.
//        OptionBuilder.withDescription("Compare input against full set of background genes, randomly select matched set of genes per permutation, and repeat [traditional bg correction]");
//        OptionBuilder.withLongOpt("properpermute");
//        OPTIONS.addOption(OptionBuilder.create("z"));
//
//        OptionBuilder.withDescription("Use logistic regression + permutation for GSEA.");
//        OptionBuilder.withLongOpt("lrt");
//        OPTIONS.addOption(OptionBuilder.create("lrt"));
//
//        OptionBuilder.withDescription("Convert input to squared values");
//        OptionBuilder.withLongOpt("square");
//        OPTIONS.addOption(OptionBuilder.create("x"));


    }

    public static void printHelp() {
//        System.out.println();
//        System.out.println("Different options for background selection exist:\n" +
//                "When using -l, foreground and background genes are limited to the set of genes specified. This option also applies to -r and -e.\n" +
//                "When using -r, a random set of background genes is selected with equal size to genes specified with -g\n" +
//                "When using -e, an equal number of background genes as specified with -g are selected by matching on expression values");
        System.out.println();
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp(" ", OPTIONS);
    }
}
