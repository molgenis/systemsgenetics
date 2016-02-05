/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.util;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.fasta.ReferenceGenomeFasta;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author MarcJan
 */
public class WriteFastaOfInterestRegion {
    
        public static void main(String[] args) {
            String referenceFasta = args[0];
            String fileWithPositions = args[1];
            String outputFile = args[2];
            
            ReferenceGenomeFasta refGen = null;
            HashMap<String, Triple<String ,Integer, Integer>> interestStrings = null;
            try {
                refGen = new ReferenceGenomeFasta(new File(referenceFasta));
            } catch (Exception ex) {
                Logger.getLogger(WriteFastaOfInterestRegion.class.getName()).log(Level.SEVERE, null, ex);
            }
            
            try {
             interestStrings= readFileWithPositions(fileWithPositions);
            } catch (Exception ex) {
                Logger.getLogger(WriteFastaOfInterestRegion.class.getName()).log(Level.SEVERE, null, ex);
            }
            
            if(refGen !=null && interestStrings !=null){
                writeFastas(refGen, interestStrings, outputFile);
            }
            
        }

    private static HashMap<String, Triple<String ,Integer, Integer>> readFileWithPositions(String fileWithPositions) throws IOException {
        HashMap<String, Triple<String ,Integer, Integer>> positions = new HashMap<>();
        TextFile gffFileReader = new TextFile(fileWithPositions, TextFile.R);
        
        String str;
        while((str=gffFileReader.readLine())!=null){
            String[] parts = str.split("\t");
            
        }
        
        gffFileReader.close();
        return positions; 
    }

    private static void writeFastas(ReferenceGenomeFasta refGen, HashMap<String, Triple<String ,Integer, Integer>> interestStrings, String outputFile) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
}
