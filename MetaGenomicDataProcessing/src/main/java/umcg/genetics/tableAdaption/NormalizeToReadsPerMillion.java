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
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author MarcJan
 */
public class NormalizeToReadsPerMillion {
    
    
    public static void main(String[] args) {
        String fileNameMetaphlanTable = "D:\\UMCG\\Projects\\nonHumanReadsRNA_Seq\\Kraken\\TB_SA_Candida_kraken.merged.noEukaryot.Filtered.tsv";
        String cladeInformationFile = "D:\\UMCG\\Projects\\nonHumanReadsRNA_Seq\\Metaphlan\\CountBased\\cladeSize.txt";
//        String fileNameMetaphlanTable = "D:\\UMCG\\Projects\\nonHumanReadsRNA_Seq\\Kraken\\TB_SA_kraken.merged.tsv";
        String fileNameReadStatistics = "D:\\UMCG\\Projects\\nonHumanReadsRNA_Seq\\ReadStatistics\\ReadStatsSAMTBLLDCandida.txt";
        
        boolean countData = true;
        boolean cladeNormalization = false;
        
        
        HashMap<String, Pair<Double, Double>> readStats = null;
        try {
            readStats = readReadStatsFile(fileNameReadStatistics);
        } catch (IOException ex) {
            Logger.getLogger(PrintTablePerLevel.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        DoubleMatrixDataset<String,String> metaphlanMatrix = null;
        try {
            metaphlanMatrix = DoubleMatrixDataset.loadDoubleData(fileNameMetaphlanTable);
        } catch (IOException ex) {
            Logger.getLogger(PrintTablePerLevel.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        if(cladeNormalization){
            HashMap<String, Double> cladeSizeInformation = null;
            try {
                cladeSizeInformation = readCladeSizeInfo(cladeInformationFile);
            } catch (IOException ex) {
                Logger.getLogger(NormalizeToReadsPerMillion.class.getName()).log(Level.SEVERE, null, ex);
            }
            metaphlanMatrix = normalizeWithCladeSize(metaphlanMatrix, cladeSizeInformation);
        }
        
        if(countData){
            metaphlanMatrix = normalizeToReadsPerMillionCount(metaphlanMatrix, readStats);
        } else {
            metaphlanMatrix = normalizeToReadsPerMillionPercent(metaphlanMatrix, readStats);
        }
        
        try {
            metaphlanMatrix.save(fileNameMetaphlanTable+".cladeNorm.normalized.tsv");
        } catch (IOException ex) {
            Logger.getLogger(PrintTablePerLevel.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    }

    private static HashMap<String, Pair<Double, Double>> readReadStatsFile(String fileNameReadStatistics) throws IOException {
        HashMap<String, Pair<Double, Double>> stats = new HashMap<String, Pair<Double, Double>>();
        
        TextFile inreader = new TextFile(fileNameReadStatistics, TextFile.R);
        
        //skip header;
        
        inreader.readLine();
        String str;
        while((str=inreader.readLine())!=null){
            String[] parts = str.split("\t");
            
            stats.put(parts[0], new Pair<Double, Double>(Double.parseDouble(parts[1]),Double.parseDouble(parts[2])));
            
//            System.out.println(str);
        }
        
        
        return stats;
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
    
    private static DoubleMatrixDataset<String, String> normalizeToReadsPerMillionPercent(DoubleMatrixDataset<String, String> metaphlanMatrix, HashMap<String, Pair<Double, Double>> readStats) {
        DoubleMatrixDataset<String, String> normalizedMatrix = new DoubleMatrixDataset<String, String>(metaphlanMatrix.getHashRows(), metaphlanMatrix.getHashCols());
        
        for(int c = 0; c<metaphlanMatrix.columns();c++){
            double normFactorPerSample = readStats.get(metaphlanMatrix.getColObjects().get(c)).getLeft();
            double normFactorPerSample2 = readStats.get(metaphlanMatrix.getColObjects().get(c)).getRight();
            for(int r = 0; r<metaphlanMatrix.rows();r++){
                normalizedMatrix.getMatrix().set(r, c, (((metaphlanMatrix.getMatrix().getQuick(r, c)*normFactorPerSample2)/normFactorPerSample)*1000000));
            }
        }
        
        return normalizedMatrix;
    }

    private static HashMap<String, Double> readCladeSizeInfo(String cladeInformationFile) throws IOException {
        HashMap<String, Double> stats = new HashMap<String, Double>();
        
        TextFile inreader = new TextFile(cladeInformationFile, TextFile.R);
        
        //skip header;
        
        inreader.readLine();
        String str;
        while((str=inreader.readLine())!=null){
            String[] parts = str.split("\t");
            
            stats.put(parts[0], Double.parseDouble(parts[1]));
            
//            System.out.println(str);
        }
        
        
        return stats;
    }

    private static DoubleMatrixDataset<String, String> normalizeWithCladeSize(DoubleMatrixDataset<String, String> metaphlanMatrix, HashMap<String, Double> cladeSizeInformation) {
        DoubleMatrixDataset<String, String> normalizedMatrix = new DoubleMatrixDataset<String, String>(metaphlanMatrix.getHashRows(), metaphlanMatrix.getHashCols());
        
        for(int r = 0; r<metaphlanMatrix.rows();r++){
            double normFactorPerClade = cladeSizeInformation.get(metaphlanMatrix.getRowObjects().get(r));
            for(int c = 0; c<metaphlanMatrix.columns();c++){
                normalizedMatrix.getMatrix().set(r, c, ((metaphlanMatrix.getMatrix().getQuick(r, c)/normFactorPerClade)));
            }
        }
        
        return normalizedMatrix;
    }

}
