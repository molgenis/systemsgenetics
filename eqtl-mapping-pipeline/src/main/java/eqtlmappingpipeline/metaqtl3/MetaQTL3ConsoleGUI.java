/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl3;

import umcg.genetica.console.ConsoleGUIElems;

/**
 *
 * @author harmjan
 */
public class MetaQTL3ConsoleGUI {

    public MetaQTL3ConsoleGUI(String[] args) {

        String settingsfile = null;
        String settingstexttoreplace = null;
        String settingstexttoreplacewith = null;
        String settingstexttoreplace2 = null;
        String settingstexttoreplace2with = null;
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
        String eqtleffectstoregressout = null;

        Double outputPlotThresold = null;
        Double maf = null;
        Double hwe = null;
        Integer nrEQTLsToOutput = null;

        String snpprobecombofile = null;
        boolean skipqqplot = false;
        boolean skipdotplot = false;
        Long rSeed = System.currentTimeMillis();

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
            } else if (arg.equals("--replacetext2")) {
                settingstexttoreplace2 = val;
            } else if (arg.equals("--replacetext2with")) {
                settingstexttoreplace2with = val;
            } else if (arg.equals("--in")) {
                in = val;
            } else if (arg.equals("--out")) {
                out = val;
            } else if (arg.equals("--text")) {
                textout = true;
            } else if (arg.equals("--binary")) {
                binout = true;
            } else if (arg.equals("--plot")) {
                try {
                    outputPlotThresold = Double.parseDouble(val);
                } catch (NumberFormatException e) {
                    System.err.println("ERROR value supplied for --plot is not a double!");
                }
            } else if (arg.equals(
                    "--inexp")) {
                inexp = val;
            } else if (arg.equals(
                    "--inexpplatform")) {
                inexpplatform = val;
            } else if (arg.equals(
                    "--inexpannot")) {
                inexpannot = val;
            } else if (arg.equals(
                    "--gte")) {
                gte = val;
            } else if (arg.equals(
                    "--cis")) {
                cis = true;
            } else if (arg.equals(
                    "--trans")) {
                trans = true;
            } else if (arg.equals(
                    "--snps")) {
                snpfile = val;
            } else if (arg.equals(
                    "--regressouteqtls")) {
                eqtleffectstoregressout = val;

            } else if (arg.equals(
                    "--snpprobe")) {
                snpprobecombofile = val;

            } else if (arg.equals(
                    "--perm")) {
                try {
                    perm = Integer.parseInt(val);
                } catch (NumberFormatException e) {
                    System.out.println("Please supply an integer for --perm");
                }
            } else if (arg.equals(
                    "--maf")) {
                try {
                    maf = Double.parseDouble(val);
                } catch (NumberFormatException e) {
                    System.out.println("Please supply an integer for --maf");
                }
            }  else if (arg.equals(
                    "--hwe")) {
                try {
                    hwe = Double.parseDouble(val);
                } catch (NumberFormatException e) {
                    System.out.println("Please supply an integer for --hwe");
                }
            } else if (arg.equals(
                    "--threads")) {
                try {
                    threads = Integer.parseInt(val);
                } catch (NumberFormatException e) {
                    System.err.println("Error --threads should be an integer");
                }

            } else if (arg.equals(
                    "--maxresults")) {
                try {
                    nrEQTLsToOutput = Integer.parseInt(val);
                } catch (NumberFormatException e) {
                    System.err.println("Error --maxresults should be an integer");
                }

            } else if (arg.equals(
                    "--skipdotplot")) {
                skipdotplot = true;
            } else if (arg.equals(
                    "--skipqqplot")) {
                skipqqplot = true;
            } else if (arg.equals(
                    "--rseed")) {
                try {
                    rSeed = Long.parseLong(val);
                } catch (NumberFormatException e) {
                    System.err.println("Error --rseed should be an integer");
                }

            }
        }
        try {
            if (settingsfile == null && in == null) {
                System.out.println("ERROR: Please supply settings file (--settings settings.xml) or --in and --out");
                printUsage();
            } else {
                MetaQTL3 m = new MetaQTL3();
                if (!binout && !textout) {
                    textout = true;
                }
                m.initialize(settingsfile, settingstexttoreplace, settingstexttoreplacewith, settingstexttoreplace2, settingstexttoreplace2with, in, inexp, inexpplatform, inexpannot, gte, out, cis, trans, perm, textout, binout, snpfile, threads, nrEQTLsToOutput, eqtleffectstoregressout, snpprobecombofile, skipdotplot, skipqqplot, rSeed, maf, hwe);
                
                if(outputPlotThresold!=null){
                    m.setOutputPlotThreshold(outputPlotThresold);
                }
                m.mapEQTLs();
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(-1);
        }
    }

    private void printUsage() {
        System.out.print("\nMetaQTL3\n" + ConsoleGUIElems.LINE);
        System.out.println("MetaQTL3 is the main eQTL mapping program of the pipeline.");
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
                + "--maf\t\t\tint\tMinimal minor allel frequency to take SNP in the analysis\n"
                + "--hwe\t\t\tint\tHardy weinberg equilibrium for SNP selection during the analysis\n"
                + "--text\t\t\t\t\tOutput results in text format\n"
                + "--binary\t\t\t\tOutput results in binary format\n"
                + "--inexp\t\t\tstring\t\tLocation of expression data\n"
                + "--inexpplatform\t\tstring\t\tGene expression platform\n"
                + "--inexpannot\t\tstring\t\tLocation of annotation file for gene expression data\n"
                + "--gte\t\t\tstring\t\tLocation of genotype to expression coupling file\n"
                + "--snps\t\t\tstring\t\tLocation of file containing SNPs to confine to\n"
                + "--threads\t\tinteger\t\tNumber of threads to calculate with. Default is number of processors.\n"
                + "--maxresults\t\tinteger\t\tNumber of results to output.\n"
                + "--regressouteqtls\tstring\t\tRegress out these eQTL effects before starting the analysis.\n"
                + "--snpprobe\t\tstring\t\tTest only combinations of SNPs and probes.");
        System.out.println("");
    }
}
