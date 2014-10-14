

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.util;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;


/**
 *
 * @author harmjan
 */
public class QTLFileSorter {

    public void run(String efile, String outfile) throws IOException {
        System.out.println("Loading: " + efile);
        QTLTextFile e = new QTLTextFile(efile, QTLTextFile.R);

        ArrayList<EQTL> eQtls = e.readList();

        e.close();
        
        System.out.println("Loaded " + eQtls.size() + " eqtls. Now Sorting.");

        Collections.sort(eQtls);

        System.out.println("Done sorting");

        if(outfile == null){
            outfile = efile + "_sorted.txt.gz";
        }
        
        QTLTextFile out = new QTLTextFile(outfile, QTLTextFile.W);
        out.write(eQtls);
        out.close();

    }
}
