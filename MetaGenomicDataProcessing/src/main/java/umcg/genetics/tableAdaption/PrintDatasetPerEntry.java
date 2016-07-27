/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetics.tableAdaption;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.MatrixHandling;

/**
 * @author MarcJan
 */
public class PrintDatasetPerEntry {

    public static void main(String[] args) {
        
//        String inputTable = "D:\\UMCG\\Projects\\MGS_MicrobiomeQTLs\\Kegg_Genes\\LLD_geneFamilies_KeggTable_KO.tsv.QuantileNormalized.txt";
//        String outputFolder = "D:\\UMCG\\Projects\\MGS_MicrobiomeQTLs\\Kegg_Genes\\PerEntry\\";

        String inputTable = "D:\\OnlineFolders\\Dropbox\\Dropbox\\miQTL\\bacteriaLogit_2.txt";
        String outputFolder = "D:\\OnlineFolders\\Dropbox\\Dropbox\\miQTL\\PerEntry\\";

//        String inputTable = "E:\\OnlineFolders\\BitSync\\UMCG\\Projects\\MGS_MicrobiomeQTLs\\GFD_IBS_IBD_LLD_GoNL_500Fg_metaphlan_2.2_results_AsinNorm2.txt";
//        String outputFolder = "E:\\OnlineFolders\\BitSync\\UMCG\\Projects\\MGS_MicrobiomeQTLs\\DataPerEntry\\";
      
//        String inputTable = "D:\\UMCG\\Projects\\MGS_MicrobiomeQTLs\\Pathways\\DataPerCohort\\GFD_LLD_500Fg_IBS_IBD_pathways_selected.tsv.QuantileNormalized.ProbesCentered.txt";
//        String outputFolder = "D:\\UMCG\\Projects\\MGS_MicrobiomeQTLs\\Pathways\\DataPerEntry\\";
        
//        String inputTable = "D:\\UMCG\\Projects\\MGS_MicrobiomeQTLs\\GO_Grouping\\AllMGS_KO50000.tsv.QuantileNormalized.txt.gz";
//        String outputFolder = "D:\\UMCG\\Projects\\MGS_MicrobiomeQTLs\\GO_Grouping\\DataPerEntry_KO50000\\";

        int cutOffNumber = 5;
        int minPercentage = 0;
        boolean percentage = false;
        boolean takeCuttOffAndPercentage = false;
        boolean removeStrainInformation = false;
        
        DoubleMatrixDataset<String, String> bugMatrix = null;
        try {
            bugMatrix = DoubleMatrixDataset.loadDoubleData(inputTable);
        } catch (IOException ex) {
            Logger.getLogger(PrintDatasetPerEntry.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        //Write matrix with direct filtering on X number of bugs in the individual cohort.
        writeTablesToFile(bugMatrix, cutOffNumber, minPercentage, percentage, takeCuttOffAndPercentage, removeStrainInformation, outputFolder);
        
    }
    
    private static void writeTablesToFile(DoubleMatrixDataset<String, String> bugMatrix, int cutOffNumber, int minPercentage, boolean percentage, boolean takeCuttOffAndPercentage, boolean removeStrainInformation, String outputFolder) {
        int i=0;
        for(String row : bugMatrix.getRowObjects()){
            // select samples of interest
            HashSet<String> s = new HashSet<String>();
            s.add(row);
            DoubleMatrixDataset<String, String> relevantBugMatrix = MatrixHandling.CreatSubsetBasedOnRows(bugMatrix, s, false);

            //Select bugs present in at least cutOffNumer
            relevantBugMatrix =  selectToppresentBugs(relevantBugMatrix,cutOffNumber, minPercentage, percentage, takeCuttOffAndPercentage, removeStrainInformation);
            
            if(relevantBugMatrix!=null){
                String[] names = row.split("\\|");
                String partToWrite = names[names.length-1];
                    
                try {
                    relevantBugMatrix.save(outputFolder+"matrix_"+i+"_CuttOfNumber_"+cutOffNumber+"_incPercentage"+minPercentage+".tsv");
                } catch (IOException ex) {
                    Logger.getLogger(PrintDatasetPerEntry.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
            i++;
        }
        
    }

    private static DoubleMatrixDataset<String, String> selectToppresentBugs(DoubleMatrixDataset<String, String> relevantBugMatrix, int cutOffNumber, int minPercentage, boolean percentage, boolean takeCuttOffAndPercentage, boolean removeStrainInformation) {
        HashSet<String> interestingColumns = new HashSet<String>();
        double nrSamples = relevantBugMatrix.columns();
        for(int r = 0; r < relevantBugMatrix.rows(); r++){
                int nonZeroCols = 0;
                for(int c = 0; c < nrSamples; c++){
                    if(relevantBugMatrix.getMatrix().getQuick(r, c)>0){
                        nonZeroCols++;
                        interestingColumns.add(relevantBugMatrix.getColObjects().get(c));
                    }
                }
                if(!percentage && nonZeroCols>=cutOffNumber){
                    return MatrixHandling.CreatSubsetBasedOnColumns(relevantBugMatrix, interestingColumns, false);
                } else if(takeCuttOffAndPercentage && ((nonZeroCols/nrSamples)*100)>=minPercentage && nonZeroCols>=cutOffNumber){
                    return MatrixHandling.CreatSubsetBasedOnColumns(relevantBugMatrix, interestingColumns, false);
                } else if(percentage && ((nonZeroCols/nrSamples)*100)>=minPercentage && nonZeroCols>=cutOffNumber){
                    return MatrixHandling.CreatSubsetBasedOnColumns(relevantBugMatrix, interestingColumns, false);
                }
        }
        return null;
        
    }
    
    public static void listFilesForFolder(File folder, HashSet<File> fileList) {
        for (final File fileEntry : folder.listFiles()) {
            if (fileEntry.isDirectory()) {
                listFilesForFolder(fileEntry, fileList);
            } else {
                fileList.add(fileEntry);
            }
        }
    }
}
