/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetics.tableAdaption;

import java.io.IOException;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.containers.Pair;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author MarcJan
 */
public class NormalizeAndFilter {
    
    
    public static void main(String[] args) {
        String fileNameMetaphlanTable = "D:\\UMCG\\Projects\\nonHumanReadsRNA_Seq\\Metaphlan\\CountBased\\Merged_metaphlan_2.2_results_LLD_TB_SA_reads.tsv";
        String additionToFileName="";
        
        boolean kraken = true;
        boolean scaleToPercent = true;
        boolean toArcSin = true;
        
//        boolean filterLowAbundantReads = true;
//        boolean filterLowPresence = true;
        
        
        DoubleMatrixDataset<String,String> metaphlanMatrix = null;
        try {
            metaphlanMatrix = DoubleMatrixDataset.loadDoubleData(fileNameMetaphlanTable);
        } catch (IOException ex) {
            Logger.getLogger(PrintTablePerLevel.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        if(scaleToPercent){
            additionToFileName+=".percentage";
            if(kraken){
                metaphlanMatrix = normalizeToPercentKraken(metaphlanMatrix);
            } else {
                metaphlanMatrix = normalizeToPercentMetaphlan(metaphlanMatrix);
            }
            
        }
        
        if(toArcSin){
            additionToFileName+=".asin";
            metaphlanMatrix = normalizeAsinSqrt(metaphlanMatrix);
        } 
        
        try {
            metaphlanMatrix.save(fileNameMetaphlanTable+additionToFileName+".tsv");
        } catch (IOException ex) {
            Logger.getLogger(PrintTablePerLevel.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    }

    private static DoubleMatrixDataset<String, String> normalizeToReadsPerMillionCount(DoubleMatrixDataset<String, String> metaphlanMatrix, HashMap<String, Pair<Double, Double>> readStats) {
        DoubleMatrixDataset<String, String> normalizedMatrix = new DoubleMatrixDataset<String, String>(metaphlanMatrix.getHashRows(), metaphlanMatrix.getHashCols());
        
        for(int c = 0; c<metaphlanMatrix.columns();c++){
            double normFactorPerSample = readStats.get(metaphlanMatrix.getColObjects().get(c)).getLeft();
            for(int r = 0; r<metaphlanMatrix.rows();r++){
                normalizedMatrix.getMatrix().set(r, c, ((metaphlanMatrix.getMatrix().getQuick(r, c)/normFactorPerSample)*1000000));
            }
        }
        
        return normalizedMatrix;
    }

    private static DoubleMatrixDataset<String, String> normalizeToPercentKraken(DoubleMatrixDataset<String, String> metaphlanMatrix) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    private static DoubleMatrixDataset<String, String> normalizeToPercentMetaphlan(DoubleMatrixDataset<String, String> metaphlanMatrix) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    private static DoubleMatrixDataset<String, String> normalizeAsinSqrt(DoubleMatrixDataset<String, String> metaphlanMatrix) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    

}
