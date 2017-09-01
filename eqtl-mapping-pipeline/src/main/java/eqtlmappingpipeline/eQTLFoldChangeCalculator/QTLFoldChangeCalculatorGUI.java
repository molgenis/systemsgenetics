/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.eQTLFoldChangeCalculator;

import umcg.genetica.console.ConsoleGUIElems;

/**
 *
 * @author harmjan
 */
public class QTLFoldChangeCalculatorGUI {

    public QTLFoldChangeCalculatorGUI(String[] args) {

        String settingsfile = null;
        String settingstexttoreplace = null;
        String settingstexttoreplacewith = null;
        String in = null;
        String out = null;
        boolean cis = false;
        boolean trans = false;
        int perm = 1;
        String outtype = "text";
        String inexp = null;
        String inexpplatform = null;
        String inexpannot = null;
        String gte = null;
        String snpfile = null;
        Integer threads = null;
        boolean textout = false;
        boolean binout = false;
        String eqtlfile = null;
        Double maf = null;
        Double hwe = null;

        for (int i = 0; i < args.length; i++) {
            String arg = args[i];
            String val = null;

            if (i + 1 < args.length) {
                val = args[i + 1];
            }

            if (arg.equals("--settings")) {
                settingsfile = val;
            } else if (arg.equals("--replacetext")) {
                settingstexttoreplace = val;
            } else if (arg.equals("--replacetextwith")) {
                settingstexttoreplacewith = val;
            } else if (arg.equals("--in")) {
                in = val;
            } else if (arg.equals("--out")) {
                out = val;
            } else if (arg.equals("--text")) {
                textout = true;
            } else if (arg.equals("--binary")) {
                binout = true;
            } else if (arg.equals("--inexp")) {
                inexp = val;
            } else if (arg.equals("--inexpplatform")) {
                inexpplatform = val;
            } else if (arg.equals("--inexpannot")) {
                inexpannot = val;
            } else if (arg.equals("--gte")) {
                gte = val;
            } else if (arg.equals("--cis")) {
                cis = true;
            } else if (arg.equals("--trans")) {
                trans = true;
            } else if (arg.equals("--snps")) {
                snpfile = val;
            } else if (arg.equals("--eqtls")) {
                eqtlfile = val;
            } else if (arg.equals(
                    "--maf")) {
                try {
                    maf = Double.parseDouble(val);
                } catch (NumberFormatException e) {
                    System.out.println("Please supply an integer for --maf");
                }
            } else if (arg.equals(
                    "--hwe")) {
                try {
                    hwe = Double.parseDouble(val);
                } catch (NumberFormatException e) {
                    System.out.println("Please supply an integer for --hwe");
                }
            } else if (arg.equals("--perm")) {
                try {
                    perm = Integer.parseInt(val);
                } catch (NumberFormatException e) {
                    System.out.println("Please supply an integer for --perm");
                }
            } else if (arg.equals("--threads")) {
                try {
                    threads = Integer.parseInt(val);
                } catch (NumberFormatException e) {
                    System.err.println("Error --threads should be an integer");
                }

            }
        }

        try {
            if (settingsfile == null && in == null) {
                System.out.println("ERROR: Please supply settings file (--settings settings.xml) or --in and --out");
                printUsage();
            } else {
                QTLFoldChangeCalculator m = new QTLFoldChangeCalculator();
                if (!binout && !textout) {
                    textout = true;
                }
                m.initialize(settingsfile, settingstexttoreplace, settingstexttoreplacewith, null, null, in, inexp, inexpplatform, inexpannot, gte, out, cis, trans, perm, textout, binout, snpfile, threads, null, null, null, true, true, null, maf, hwe);
                m.calculateFoldChanges(eqtlfile);
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(-1);
        }
    }

    private void printUsage() {
        System.out.print("\nFoldchange\n" + ConsoleGUIElems.LINE);
        System.out.println("Foldchange calculation for a list of eQTL effects (statistically not very useful, but included for 'legacy' purposes.");
        System.out.print("\nExamples\n" + ConsoleGUIElems.LINE);
        System.out.println("Example using settingsfile:\tjava -jar eQTLMappingPipeline.jar --mode metaqtl --settings settings.xml");
        System.out.println("Example using commandline:\tjava -jar eQTLMappingPipeline.jar --mode metaqtl --in /path/to/GenotypeMatrix.dat --out /path/to/output/ --cis --perm 10 --text --inexp /path/to/expressiondata.txt --inexpannot /path/to/annotation.txt --gte /path/to/genotypetoexpressioncoupling.txt");
        System.out.println("");
        System.out.print("Settings file options:\n" + ConsoleGUIElems.LINE);
        System.out.println("--settings\t\tsettings.xml\tLocation of settings file\n"
                + "--replacetext\t\ttext\t\tText to replace in settings file\n"
                + "--replacetextwith\ttext\t\tReplace the text in the settings file, defined by --replacetext with the following text (can be empty)");

        System.out.println("");
        System.out.print("Command line options:\n" + ConsoleGUIElems.LINE);
        System.out.println("--in\t\t\tdir\t\tLocation of the genotype data\n"
                + "--out\t\t\tdir\t\tLocation where the output should be stored\n"
                + "--cis\t\t\t\t\tPerform cis-eQTL analysis\n"
                + "--trans\t\t\t\t\tPerform trans-eQTL analysis\n"
                + "--perm\t\t\tint\t\tNumber of permutations to perform\n"
                + "--text\t\t\t\t\tOutput results in text format\n"
                + "--binary\t\t\t\tOutput results in binary format\n"
                + "--inexp\t\t\tstring\t\tLocation of expression data\n"
                + "--inexpplatform\t\tstring\t\tGene expression platform\n"
                + "--inexpannot\t\tstring\t\tLocation of annotation file for gene expression data\n"
                + "--gte\t\t\tstring\t\tLocation of genotype to expression coupling file\n"
                + "--snps\t\t\tstring\t\tLocation of file containing SNPs to confine to\n"
                + "--threads\t\tinteger\t\tNumber of threads to calculate with. Default is number of processors."
                + "--eqtls\t\tFile containing the eQTLs");
        System.out.println("");
    }
}
