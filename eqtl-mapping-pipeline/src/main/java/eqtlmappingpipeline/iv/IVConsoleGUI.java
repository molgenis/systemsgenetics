/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.iv;

import eqtlmappingpipeline.metaqtl3.MetaQTL3;
import umcg.genetica.console.ConsoleGUIElems;

/**
 *
 * @author harmjan
 */
public class IVConsoleGUI {

    public IVConsoleGUI(String[] args) {
        String settingsfile = null;
        
        String in = null;
        String out = null;
        boolean cis = false;
        boolean trans = false;
        int perm = 1;
        
        String inexp = null;
        String inexpplatform = null;
        String inexpannot = null;
        String gte = null;
        String snpProbeCombinationList = null;

        for (int i = 0; i < args.length; i++) {
            String arg = args[i];
            String val = null;

            if (i + 1 < args.length) {
                val = args[i + 1];
            }

            if (arg.equals("--settings")) {
                settingsfile = val;
            } else if (arg.equals("--in")) {
                in = val;
            } else if (arg.equals("--out")) {
                out = val;
            } else if (arg.equals("--inexp")) {
                inexp = val;
            } else if (arg.equals("--inexpplatform")) {
                inexpplatform = val;
            } else if (arg.equals("--inexpannot")) {
                inexpannot = val;
            } else if (arg.equals("--gte")) {
                gte = val;
            } else if (arg.equals("--effects")) {
                snpProbeCombinationList = val;
            } else if (arg.equals("--perm")) {
                try {
                    perm = Integer.parseInt(val);
                } catch (NumberFormatException e) {
                    System.out.println("Please supply an integer for --perm");
                }
            } 
        }

        try {
            if (settingsfile == null && in == null) {
                System.out.println("ERROR: Please supply settings file (--settings settings.xml) or --in and --out");
                printUsage();
            } else {
                new IVAnalysis(settingsfile, in, inexp, inexpplatform, inexpannot, gte, out, perm, snpProbeCombinationList);
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(-1);
        }
    }

    private void printUsage() {
        System.out.print("\nIVAnalysis\n" + ConsoleGUIElems.LINE);
        System.out.println("IVAnalysis uses instrumental variables in regression, giving a gain in power towards trans-eQTL mapping and possible clues on causality");
        System.out.print("\nExamples\n" + ConsoleGUIElems.LINE);
        System.out.println("Example using settingsfile:\tjava -jar eQTLMappingPipeline.jar --mode iv --settings settings.xml");
        System.out.println("Example using commandline:\tjava -jar eQTLMappingPipeline.jar --mode iv --in /path/to/GenotypeMatrix.dat --out /path/to/output/ --perm 10 --inexp /path/to/expressiondata.txt --inexpannot /path/to/annotation.txt --gte /path/to/genotypetoexpressioncoupling.txt");
        System.out.println("");
        System.out.print("Settings file options:\n" + ConsoleGUIElems.LINE);
        System.out.println("--settings\t\tsettings.xml\tLocation of settings file\n");

        System.out.println("");
        System.out.print("Command line options:\n" + ConsoleGUIElems.LINE);
        System.out.println("--in\t\t\tdir\t\tLocation of the genotype data\n"
                + "--out\t\t\tdir\t\tLocation where the output should be stored\n"
                + "--perm\t\t\tint\t\tNumber of permutations to perform\n"
                + "--inexp\t\t\tstring\t\tLocation of expression data\n"
                + "--inexpplatform\t\tstring\t\tGene expression platform\n"
                + "--inexpannot\t\tstring\t\tLocation of annotation file for gene expression data\n"
                + "--gte\t\t\tstring\t\tLocation of genotype to expression coupling file\n"
                + "--effects\t\t\tstring\t\tLocation of file containing cis and trans effects per SNP. Probe IDs should be similar to ids used in expressiondata.txt. Tab separated format:  snp cis trans\n"
                );
        System.out.println("");
    }
}
