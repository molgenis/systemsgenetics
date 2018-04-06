/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetics.tableAdaption;

import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author MarcJan
 */
public class PrintTablePerLevel {

    public static void main(String[] args) {
//        String fileNameMetaphlanTable = "D:\\UMCG\\Projects\\nonHumanReadsRNA_Seq\\Kraken\\TB_SA_Candida_kraken.merged.noEukaryot.Filtered.normalized.tsv";
//        String fileNameMetaphlanTable = "D:\\UMCG\\Projects\\nonHumanReadsRNA_Seq\\Metaphlan\\CountBased\\Merged_metaphlan_2.2_results_LLD_TB_SA_Candida_reads_cleaned.cladeNorm.normalized.tsv";
//        String fileNameMetaphlanTable = "D:\\UMCG\\Projects\\nonHumanReadsRNA_Seq\\Metaphlan\\CountBased\\Merged_metaphlan_2.2_results_LLD_TB_SA_Candida_reads_cleaned.normalized.tsv";
//        String fileNameMetaphlanTable = "D:\\UMCG\\Projects\\MGS_MicrobiomeQTLs\\phenotypes\\LLD_metaphlan_2.2_results_AsinNorm.ProbesWithZeroVarianceRemoved.CovariatesRemoved.txt";
        String fileNameMetaphlanTable = "D:\\UMCG\\Projects\\MGS_MicrobiomeQTLs\\GFD_IBS_IBD_LLD_GoNL_500Fg_metaphlan_2.2_results2.txt";
        boolean kraken = false;

        DoubleMatrixDataset<String, String> metaphlanMatrix = null;
        try {
            metaphlanMatrix = DoubleMatrixDataset.loadDoubleData(fileNameMetaphlanTable);
        } catch (IOException ex) {
            Logger.getLogger(PrintTablePerLevel.class.getName()).log(Level.SEVERE, null, ex);
        }

        if (kraken) {
            try {
                writeMatrixPerLevelKraken(metaphlanMatrix, fileNameMetaphlanTable);
            } catch (IOException ex) {
                Logger.getLogger(PrintTablePerLevel.class.getName()).log(Level.SEVERE, null, ex);
            }
        } else {
            try {
                writeMatrixPerLevelMetaphlan(metaphlanMatrix, fileNameMetaphlanTable);
            } catch (IOException ex) {
                Logger.getLogger(PrintTablePerLevel.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

    }

    private static void writeMatrixPerLevelKraken(DoubleMatrixDataset<String, String> metaphlanMatrix, String fileNameMetaphlanTable) throws IOException {
        String[] levels = {"d__", "|k__", "|p__", "|c__", "|o__", "|f__", "|g__", "|s__"};
        String[] names = {"_domain", "_kingdom", "_phylum", "_class", "_order", "_family", "_genus", "_species"};

        ArrayList<String> alreadyPrintedData = new ArrayList<String>();
        for (int i = levels.length; i > 0; i--) {
//            System.out.println(levels[i-1]);
//            System.out.println(names[i-1]);
            ArrayList<String> relevantEntries = new ArrayList<String>();

            for (String rowObj : metaphlanMatrix.getRowObjects()) {
//                System.out.println(rowObj);
                if (rowObj.contains(levels[i - 1]) && !alreadyPrintedData.contains(rowObj)) {
//                    System.out.println("\t"+rowObj);
                    alreadyPrintedData.add(rowObj);
                    relevantEntries.add(rowObj);
                }
            }
            if (relevantEntries.size() != 0) {
                LinkedHashMap<String, Integer> levelEntries = new LinkedHashMap<String, Integer>();
                int rT = 0;
                for (String key : relevantEntries) {
                    levelEntries.put(key, rT);
                    rT++;
                }

                DoubleMatrixDataset<String, String> metaphlanLevelMatrix = new DoubleMatrixDataset<String, String>(levelEntries, metaphlanMatrix.getHashCols());

                for (int r = 0; r < metaphlanLevelMatrix.rows(); r++) {
                    int originalR = metaphlanMatrix.getHashRows().get(metaphlanLevelMatrix.getRowObjects().get(r));
                    for (int c = 0; c < metaphlanLevelMatrix.columns(); c++) {
                        metaphlanLevelMatrix.getMatrix().setQuick(r, c, metaphlanMatrix.getMatrix().getQuick(originalR, c));
                    }
                }

                metaphlanLevelMatrix.save(fileNameMetaphlanTable + names[i - 1] + ".tsv");
            }
        }
    }

    private static void writeMatrixPerLevelMetaphlan(DoubleMatrixDataset<String, String> metaphlanMatrix, String fileNameMetaphlanTable) throws IOException {
        String[] levels = {"k__", "|p__", "|c__", "|o__", "|f__", "|g__", "|s__", "|t__", "_unclassified"};
        String[] names = {"_kingdom", "_phylum", "_class", "_order", "_family", "_genus", "_species", "_strain", "_endPoints"};

        ArrayList<String> alreadyPrintedData = new ArrayList<String>();
        for (int i = levels.length; i > 0; i--) {
//            System.out.println(levels[i-1]);
//            System.out.println(names[i-1]);
            ArrayList<String> relevantEntries = new ArrayList<String>();
            if (i == levels.length) {
                for (String rowObj : metaphlanMatrix.getRowObjects()) {
                    //                System.out.println(rowObj);
                    if (rowObj.endsWith(levels[i - 1]) || rowObj.contains(levels[i - 2])) {
//                        System.out.println("\t"+rowObj);
                        relevantEntries.add(rowObj);
                    }
                }
            } else {
                for (String rowObj : metaphlanMatrix.getRowObjects()) {
                    //                System.out.println(rowObj);
                    if (rowObj.contains(levels[i - 1]) && !alreadyPrintedData.contains(rowObj)) {
//                        System.out.println("\t"+rowObj);
                        alreadyPrintedData.add(rowObj);
                        relevantEntries.add(rowObj);
                    }
                }
            }

            if (relevantEntries.size() != 0) {
                LinkedHashMap<String, Integer> levelEntries = new LinkedHashMap<String, Integer>();
                int rT = 0;
                for (String key : relevantEntries) {
                    levelEntries.put(key, rT);
                    rT++;
                }

                DoubleMatrixDataset<String, String> metaphlanLevelMatrix = new DoubleMatrixDataset<String, String>(levelEntries, metaphlanMatrix.getHashCols());

                for (int r = 0; r < metaphlanLevelMatrix.rows(); r++) {
                    int originalR = metaphlanMatrix.getHashRows().get(metaphlanLevelMatrix.getRowObjects().get(r));
                    for (int c = 0; c < metaphlanLevelMatrix.columns(); c++) {
                        metaphlanLevelMatrix.getMatrix().setQuick(r, c, metaphlanMatrix.getMatrix().getQuick(originalR, c));
                    }
                }

                metaphlanLevelMatrix.save(fileNameMetaphlanTable + names[i - 1] + ".tsv");
            }
        }
    }
}
