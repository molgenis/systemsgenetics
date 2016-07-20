/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetics.tableAdaption;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author MarcJan
 */
public class ExtendSnpProbeCombinationsFile {
    
    public static void main(String[] args) {
        String fileNameMetaphlanTable = "D:\\UMCG\\Projects\\MGS_MicrobiomeQTLs\\GFD_IBS_IBD_LLD_GoNL_500Fg_metaphlan_2.2_results2.txt";
        String initialReplicationList = "D:\\OnlineFolders\\AeroFS\\SharedQTLFolder\\QTL_Output\\For_paper_20160208\\LifeLines_Discovery\\bQTL_Replist_nonZero.txt";
        String fullReplicationList = "D:\\OnlineFolders\\AeroFS\\SharedQTLFolder\\QTL_Output\\For_paper_20160208\\LifeLines_Discovery\\bQTL_Replist_Full_nonZero.txt";
        
        DoubleMatrixDataset<String, String> metaphlanMatrix = null;
        ArrayList<Pair<String,String>> initialPairs = null;
        
        try {
            initialPairs = readInitialPairs(initialReplicationList);
        } catch (IOException ex) {
            Logger.getLogger(ExtendSnpProbeCombinationsFile.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        try {
            metaphlanMatrix = DoubleMatrixDataset.loadDoubleData(fileNameMetaphlanTable);
        } catch (IOException ex) {
            Logger.getLogger(PrintTablePerLevel.class.getName()).log(Level.SEVERE, null, ex);
        }
        if(metaphlanMatrix!=null && initialPairs!=null){
            ArrayList<String> rows = metaphlanMatrix.getRowObjects();

            try {
                writeFullReplicationList(rows, initialPairs, fullReplicationList);
            } catch (IOException ex) {
                Logger.getLogger(ExtendSnpProbeCombinationsFile.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

    }

    private static ArrayList<Pair<String, String>> readInitialPairs(String initialReplicationList) throws IOException {
        ArrayList<Pair<String, String>> combinationsToTest = new ArrayList<Pair<String, String>>();
        TextFile reader = new TextFile(initialReplicationList, TextFile.R);
        
        String st;
        while((st=reader.readLine())!=null){
            String[] parts = st.split("\t");
            combinationsToTest.add(new Pair<String, String>(parts[0], parts[1]));
        }
        
        reader.close();
        return combinationsToTest;
    }

    private static void writeFullReplicationList(ArrayList<String> rows, ArrayList<Pair<String, String>> initialPairs, String fullReplicationList) throws IOException {
        TextFile writer = new TextFile(fullReplicationList, TextFile.W);
        
        HashSet<String> previousPairs = new HashSet<String>();
        for(Pair<String, String> p : initialPairs){
            for(String b: rows){
                StringBuilder combi = new StringBuilder();
                combi.append(p.getLeft()).append('\t').append(b);
                if(b.contains(p.getRight()) && !previousPairs.contains(combi.toString())){
                    writer.writeln(combi.toString());
                    previousPairs.add(combi.toString());
                }
            }
        }
        
        writer.close();
    }
}
