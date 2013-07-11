/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.util.eqtlfilesorter;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class EQTLFileSorter {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        EQTLFileSorter ef = new EQTLFileSorter();
        try {
            for (int perm = 0; perm < 11; perm++) {
                String dir = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQCReWeightedZ-Batch2/";
                String sortedDir = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQCReWeightedZ-Batch2Sorted/";


                String infile = "eQTLs.txt.gz";
                if (perm > 0) {
                    infile = "PermutedEQTLsPermutationRound" + perm + ".txt.gz";

                }
                ef.run(dir + infile, null);

                Gpio.createDir(sortedDir);

                Gpio.moveFile(dir + infile + "_sorted.txt.gz", sortedDir + infile);
            }

//	    String filename = "/Volumes/iSnackHD/Data/Projects/Isis/eQTLResults/2012-10-30-eQTLs-Cis/SNP-Initial/QC/eQTLs.txt";
//	    ef.run(filename);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void run(String efile, String outfile) throws IOException {
        System.out.println("Loading: " + efile);
        TextFile e = new TextFile(efile, TextFile.R);

        String[] elems = e.readLineElemsReturnObjects(TextFile.tab);

        String header = Strings.concat(elems, Strings.tab);

        ArrayList<SortableEQTL> eqtls = new ArrayList<SortableEQTL>();

        elems = e.readLineElemsReturnObjects(TextFile.tab);
        int lnctr = 0;
        while (elems != null) {

            SortableEQTL eqtl = new SortableEQTL();
            eqtl.absZScore = Math.abs(Double.parseDouble(elems[10]));
            eqtl.line = Strings.concat(elems, Strings.tab);
            eqtl.pvalue = Double.parseDouble(elems[0]);
            eqtls.add(eqtl);
            elems = e.readLineElemsReturnObjects(TextFile.tab);
            lnctr++;
            if (lnctr % 1000000 == 0) {
                System.out.println(lnctr + " lines parsed");

            }
        }



        e.close();
        System.out.println("Loaded " + eqtls.size() + " eqtls. Now Sorting.");

        Collections.sort(eqtls);

        System.out.println("Done sorting");

        if(outfile == null){
            outfile = efile + "_sorted.txt.gz";
        }
        
        TextFile out = new TextFile(outfile, TextFile.W);
        out.writeln(header);
        out.writeList(eqtls);
        out.close();

//	FDR.generateEProbesFile(efile + "_sorted.txt", efile + "eQTLProbesFDR0.05_sorted.txt");
//	FDR.generateESNPsFile(efile + "_sorted.txt", efile + "eQTLSNPsFDR0.05_sorted.txt");
    }
}
