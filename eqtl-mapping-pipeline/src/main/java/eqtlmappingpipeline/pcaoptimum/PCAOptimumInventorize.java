/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.pcaoptimum;

import eqtlmappingpipeline.util.eQTLFileCompare;
import java.io.IOException;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.eQTLTextFile;

/**
 *
 * @author harmjan
 */
public class PCAOptimumInventorize {

    public void inventory(String in, boolean cis, boolean trans, int max, int stepSize) throws IOException, Exception {
        int[] nrCISEQTLsPerRound = new int[21];
        int[] nrCisSharedPerRound = new int[21];
        int[] nrCisOppositePerRound = new int[21];

        int[] nrTransEQTLsPerRound = new int[21];
        int[] nrTransSharedPerRound = new int[21];
        int[] nrTransOppositePerRound = new int[21];

        eQTLFileCompare e = new eQTLFileCompare();
        int round = 1;
        
        if(cis){
            eQTLTextFile etf = new eQTLTextFile(in + "Cis-0PCAsRemoved/eQTLProbesFDR0.05.txt", eQTLTextFile.R);
            EQTL[] origEQTLs = etf.read();
            nrCISEQTLsPerRound[0] = origEQTLs.length;
            etf.close();
        }
        if(trans){
            eQTLTextFile etf = new eQTLTextFile(in + "Trans-0PCAsRemoved/eQTLProbesFDR0.05.txt", eQTLTextFile.R);
            EQTL[] origEQTLs = etf.read();
            nrTransEQTLsPerRound[0] = origEQTLs.length;
            etf.close();
        }
        for (int pca = stepSize; pca <= max; pca += stepSize) {

            if(cis){
                String cisNull = in + "Cis-0PCAsRemoved/eQTLsFDR0.05.txt";
                String cisOut = in + "Cis-" + pca + "PCAsRemoved/eQTLsFDR0.05.txt";

                e.compareOverlapAndZScoreDirectionTwoEQTLFiles(cisOut, cisNull, in + "Cis-" + pca + "PCAsRemoved", false);
                nrCisSharedPerRound[round] = e.getNrShared();
                nrCisOppositePerRound[round] = e.getNrOpposite();

                eQTLTextFile etf2 = new eQTLTextFile(in + "Cis-" + pca + "PCAsRemoved/eQTLProbesFDR0.05.txt", eQTLTextFile.R);
                nrCISEQTLsPerRound[round] = etf2.read().length;
                etf2.close();
            }  
            
            if(trans){
                String transNull = in + "Trans-0PCAsRemoved/eQTLsFDR0.05.txt";
                String transOut = in + "Trans-" + pca + "PCAsRemoved/eQTLsFDR0.05.txt";

                e.compareOverlapAndZScoreDirectionTwoEQTLFiles(transOut, transNull, in + "Trans-" + pca + "PCAsRemoved", false);
                nrTransSharedPerRound[round] = e.getNrShared();
                nrTransOppositePerRound[round] = e.getNrOpposite();

                eQTLTextFile etf2 = new eQTLTextFile(in + "Trans-" + pca + "PCAsRemoved/eQTLProbesFDR0.05.txt", eQTLTextFile.R);
                nrTransEQTLsPerRound[round] = etf2.read().length;
                etf2.close();
            }
            round++;
        }

        System.out.println("PCs\tCis\tShared\tDifferentAllelicDirection\tTrans\tShared\tDifferentAllelicDirection");
        for (int i = 1; i <= (max/stepSize); i++) {
            System.out.println((i * stepSize) + "\t" + nrCISEQTLsPerRound[i] + "\t" + nrCisSharedPerRound[i] + "\t" + nrCisOppositePerRound[i] + "\t" + nrTransEQTLsPerRound[i] + "\t" + nrTransSharedPerRound[i] + "\t" + nrTransOppositePerRound[i]);
        }
    }

    void inventorypcqtl(String in, boolean cis, boolean trans, int max, int stepSize) throws IOException, Exception {
        int[] nrCISEQTLsPerRound = new int[21];
        int[] nrCisSharedPerRound = new int[21];
        int[] nrCisOppositePerRound = new int[21];

        int[] nrTransEQTLsPerRound = new int[21];
        int[] nrTransSharedPerRound = new int[21];
        int[] nrTransOppositePerRound = new int[21];
        eQTLFileCompare e = new eQTLFileCompare();
        int round = 1;

        String pcqtlsuffix = "-GeneticVectorsNotRemoved";
        if(cis){
            eQTLTextFile etf = new eQTLTextFile(in + "Cis-0PCAsRemoved/eQTLProbesFDR0.05.txt", eQTLTextFile.R);
            EQTL[] origEQTLs = etf.read();
            nrCISEQTLsPerRound[0] = origEQTLs.length;
            etf.close();
        }
        if(trans){
            eQTLTextFile etf = new eQTLTextFile(in + "Trans-0PCAsRemoved/eQTLProbesFDR0.05.txt", eQTLTextFile.R);
            EQTL[] origEQTLs = etf.read();
            nrTransEQTLsPerRound[0] = origEQTLs.length;
            etf.close();
        }
        for (int pca = stepSize; pca <= max; pca += stepSize) {
            if(cis){
                String cisNull = in + "Cis-0PCAsRemoved/eQTLsFDR0.05.txt";
                String cisOut = in + "Cis-" + pca + "PCAsRemoved" + pcqtlsuffix + "/eQTLsFDR0.05.txt";

                e.compareOverlapAndZScoreDirectionTwoEQTLFiles(cisOut, cisNull, in + "Cis-" + pca + "PCAsRemoved" + pcqtlsuffix, false);
                nrCisSharedPerRound[round] = e.getNrShared();
                nrCisOppositePerRound[round] = e.getNrOpposite();

                eQTLTextFile etf2 = new eQTLTextFile(in + "Cis-" + pca + "PCAsRemoved" + pcqtlsuffix + "/eQTLProbesFDR0.05.txt", eQTLTextFile.R);
                nrCISEQTLsPerRound[round] = etf2.read().length;
                etf2.close();
            }
            if(trans){
                String transNull = in + "Trans-0PCAsRemoved/eQTLsFDR0.05.txt";
                String transOut = in + "Trans-" + pca + "PCAsRemoved" + pcqtlsuffix + "/eQTLsFDR0.05.txt";

                e.compareOverlapAndZScoreDirectionTwoEQTLFiles(transOut, transNull, in + "Trans-" + pca + "PCAsRemoved" + pcqtlsuffix, false);
                nrTransSharedPerRound[round] = e.getNrShared();
                nrTransOppositePerRound[round] = e.getNrOpposite();

                eQTLTextFile etf2 = new eQTLTextFile(in + "Trans-" + pca + "PCAsRemoved" + pcqtlsuffix + "/eQTLProbesFDR0.05.txt", eQTLTextFile.R);
                nrTransEQTLsPerRound[round] = etf2.read().length;
                etf2.close();
            }
            round++;
        }

        System.out.println("PCs\tCis\tShared\tDifferentAllelicDirection\tTrans\tShared\tDifferentAllelicDirection");
        
        for (int i = 1; i <= (max/stepSize); i++) {
            System.out.println((i * stepSize) + "\t" + nrCISEQTLsPerRound[i] + "\t" + nrCisSharedPerRound[i] + "\t" + nrCisOppositePerRound[i] + "\t" + nrTransEQTLsPerRound[i] + "\t" + nrTransSharedPerRound[i] + "\t" + nrTransOppositePerRound[i]);
//            System.out.println((i * 5) + "\t" + nrCISEQTLsPerRound[i] + "\t" + nrTransEQTLsPerRound[i]);
        }
    }
}
