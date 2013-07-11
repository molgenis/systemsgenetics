/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.trityperconverter;

import java.io.IOException;
import umcg.genetica.console.ConsoleGUIElems;
import umcg.genetica.io.trityper.converters.TriTyperToMachImputedTransposed;

/**
 *
 * @author harmjan
 */
public class ReverseTriTyperConverter {

    public ReverseTriTyperConverter(String[] args) {

        String in = null;
        String out = null;
        String snps = null;

        
        for (int i = 0; i < args.length; i++) {
            String arg = args[i];
            String val = null;

            if (i + 1 < args.length) {
                val = args[i + 1];
            }

            if (arg.equals("--in")) {
                in = val;
            } else if (arg.equals("--out")) {
                out = val;
            } else if (arg.equals("--snps")) {
                snps = val;
            }
        }

        if (in == null || out == null) {
            printUsage();
        } else {
            try {
                TriTyperToMachImputedTransposed.convert(in, out, snps);
            } catch (IOException e) {
                e.printStackTrace();

            }
        }

    }

    private void printUsage() {

        System.out.print("\nTriTyper Reverse Converter\n" + ConsoleGUIElems.LINE);
        System.out.println("Use TriTyper Reverse Converter to convert your imputed TriTyper data into text format.");
        System.out.print("\nExamples\n" + ConsoleGUIElems.LINE);
        System.out.println("Example using commandline:\tjava -jar eQTLMappingPipeline.jar --mode revconvertgt --in /path/to/trityperdata/ --out /path/to/output/ [--snps /path/to/snps.txt]");
        System.out.println("");
        System.out.print("Command line options:\n" + ConsoleGUIElems.LINE);
        System.out.println("--in\t\tdir\tLocation of (imputed) TriTyper genotype data");
        System.out.println("--out\t\tdir\tLocation to save converted data");
        System.out.println("--snps\t\tstring\tFull path to SNP file");
//        System.out.println("[--type] [pedmap|mach|machtransposed|beagle|impute]\tLocation to save converted data\n");
        System.out.println("");

    }
}
