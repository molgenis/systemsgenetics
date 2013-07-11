/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.util;

import eqtlmappingpipeline.pcaoptimum.eQTLFileCompare;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.console.ConsoleGUIElems;

/**
 *
 * @author harmjan
 */
public class eQTLResultCompare {

    /**
     * @param args the command line arguments
     */
    public eQTLResultCompare(String[] args) {


        String out = null;
        String file1 = null;
        String file2 = null;
        boolean matchOnGeneName = false;

        for (int i = 0; i < args.length; i++) {
            String arg = args[i];
            String val = null;

            if (i + 1 < args.length) {
                val = args[i + 1];
            }

            if (arg.equals("--out")) {
                out = val;
            } else if (arg.equals("--file1")) {
                file1 = val;
            } else if (arg.equals("--file2")) {
                file2 = val;
            } else if (arg.equals("--genebased")) {
                matchOnGeneName = true;
            }
        }

        if (out != null && file1 != null && file2 != null) {
            try {
                eQTLFileCompare eqfc = new eQTLFileCompare();
                eqfc.compareOverlapAndZScoreDirectionTwoEQTLFiles(file1, file2, out, matchOnGeneName);
            } catch (IOException ex) {
                Logger.getLogger(eQTLResultCompare.class.getName()).log(Level.SEVERE, null, ex);
            } catch (Exception ex) {
                Logger.getLogger(eQTLResultCompare.class.getName()).log(Level.SEVERE, null, ex);
            }
        } else {
            printUsage();
        }
    }

    private void printUsage() {
       System.out.print("QTL File comparison\n" + ConsoleGUIElems.LINE);
        System.out.println("Compares two eQTL files with each other.");

        System.out.print("Command line options:\n" + ConsoleGUIElems.LINE);

        System.out.println("--out\t\tstring\t\tOutput file name\n"
                + "--file1\t\tstring\t\tLocation of file 1\n"
                + "--file2\t\tstring\t\tLocation of file 2\n"
                + "--genebased\t\t\tPerform comparison on the basis of gene names (optional, defaults to probe based comparison)\n");
    }
}
