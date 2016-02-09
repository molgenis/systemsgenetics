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
//            String referenceFasta = args[0];
//            String fileWithPositions = args[1];
//            String outputFile = args[2];

            String referenceFasta = "D:\\UMCG\\ProbeMapping\\human_g1k_v37.fasta";
            String fileWithPositions = "D:\\UMCG\\ProbeMapping\\Info\\V71_2\\GFT\\meta-exons_v71_cut_sorted_22-05-14.gtf";
            String outputFile = "D:\\UMCG\\ProbeMapping\\Test.txt";
            
            ReferenceGenomeFasta refGen = null;
            HashMap<String, HashMap<String, Triple<String ,Integer, Integer>>> interestStrings = null;
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

    private static HashMap<String, HashMap<String, Triple<String ,Integer, Integer>>> readFileWithPositions(String fileWithPositions) throws IOException {
        HashMap<String, HashMap<String, Triple<String ,Integer, Integer>>> chromosomalProbeInformation = new HashMap<>();
        TextFile gffFileReader = new TextFile(fileWithPositions, TextFile.R);
        
        int counter = 0;
        String str;
        while((str=gffFileReader.readLine())!=null){
            String[] parts = str.split("\t");
//            System.out.println(str);
            
            if(!chromosomalProbeInformation.containsKey(parts[0])){
                chromosomalProbeInformation.put(parts[0], new HashMap<String, Triple<String ,Integer, Integer>>());
            }
            
            chromosomalProbeInformation.get(parts[0]).put(parts[8], new Triple<String ,Integer, Integer>(parts[6], Integer.parseInt(parts[3]), Integer.parseInt(parts[4])));
            counter++;
        }
        System.out.println("Number of entries read in: "+counter);
        gffFileReader.close();
        return chromosomalProbeInformation; 
    }

    private static void writeFastas(ReferenceGenomeFasta refGen, HashMap<String, HashMap<String, Triple<String ,Integer, Integer>>> interestStrings, String outputFile) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
}
