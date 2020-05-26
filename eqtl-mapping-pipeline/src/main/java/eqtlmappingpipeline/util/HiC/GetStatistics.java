/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.util.HiC;


import java.io.IOException;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author MarcJan
 */
public class GetStatistics {
    
    public static void main(String[] args){
//        String inputFileReal= "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Trans_Pc22c_CisMeQTLc_meQTLs\\RegressedOut_CisEffects_New\\HiC_Annot\\New_HiC_1Kb_E30\\eQTLsFDR0.05-ProbeLevel_BsFiltered&Filtered_HiC_LD_1kb_E30_annotated.txt";
//        String inputFilePermutation= "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Trans_Pc22c_CisMeQTLc_meQTLs\\RegressedOut_CisEffects_New\\HiC_Annot\\New_HiC_1Kb_E30\\eQTLsFDR0.05-ProbeLevel_BsFiltered&Filtered_AlternativePerm_HiC_LD_1kb_E30_annotated.txt";
//        String outputFile = "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Trans_Pc22c_CisMeQTLc_meQTLs\\RegressedOut_CisEffects_New\\HiC_Annot\\New_HiC_1Kb_E30\\eQTLsFDR0.05-ProbeLevel_BsFiltered&Filtered_HiC_LD_1kb_E30_annotated_";
        String inputFileReal= "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Trans_Pc22c_CisMeQTLc_meQTLs\\RegressedOut_CisEffects_New\\HiC_Annot\\New_HiC_5Kb_0\\eQTLsFDR0.05-ProbeLevel_BsFiltered&Filtered_HiC_LD_5kb_0_annotated.txt";
        String inputFilePermutation= "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Trans_Pc22c_CisMeQTLc_meQTLs\\RegressedOut_CisEffects_New\\HiC_Annot\\New_HiC_5Kb_0\\eQTLsFDR0.05-ProbeLevel_BsFiltered&Filtered_AlternativePerm_HiC_LD_5kb_0_annotated.txt";
        String outputFile = "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Trans_Pc22c_CisMeQTLc_meQTLs\\RegressedOut_CisEffects_New\\HiC_Annot\\New_HiC_5Kb_0\\eQTLsFDR0.05-ProbeLevel_BsFiltered&Filtered_HiC_LD_5kb_0_annotated_";
        
        Pair<HashMap<String, Pair<Integer, Integer>>, HashMap<String, Pair<Integer, Integer>>> informationPerEntry = null;
        try {
            informationPerEntry = readAndProcessPermutation(inputFilePermutation);
        } catch (IOException ex) {
            Logger.getLogger(GetStatistics.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        try {
//            processRealData(inputFileReal, informationPerEntry, outputFile);
            processRealData2(readAndProcessPermutation(inputFileReal), informationPerEntry, outputFile);
        } catch (IOException ex) {
            Logger.getLogger(GetStatistics.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    private static Pair<HashMap<String, Pair<Integer, Integer>>, HashMap<String, Pair<Integer, Integer>>> readAndProcessPermutation(String inputFilePermutation) throws IOException {
        TextFile input = new TextFile(inputFilePermutation, TextFile.R);
        HashMap<String, Pair<Integer, Integer>> cpgInfo = new HashMap<>();
        HashMap<String, Pair<Integer, Integer>> snpInfo = new HashMap<>();
        
        String row;
        while((row=input.readLine())!=null){
            String[] parts = row.split("\t");
            String[] parts2 = parts[0].split("-");
            
            if(!cpgInfo.containsKey(parts2[0])){
                if(parts[1].equals("true")){
                    cpgInfo.put(parts2[0], new Pair<>(1,1));
                } else {
                    cpgInfo.put(parts2[0], new Pair<>(1,0));
                }
            } else {
                if(parts[1].equals("true")){
                    cpgInfo.put(parts2[0], new Pair<>((cpgInfo.get(parts2[0]).getLeft()+1),(cpgInfo.get(parts2[0]).getRight()+1)));
                } else {
                    cpgInfo.put(parts2[0], new Pair<>((cpgInfo.get(parts2[0]).getLeft()+1),(cpgInfo.get(parts2[0]).getRight())));
                }
            }
            if(!snpInfo.containsKey(parts2[1])){
                if(parts[1].equals("true")){
                    snpInfo.put(parts2[1], new Pair<>(1,1));
                } else {
                    snpInfo.put(parts2[1], new Pair<>(1,0));
                }
            } else {
                if(parts[1].equals("true")){
                    snpInfo.put(parts2[1], new Pair<>((snpInfo.get(parts2[1]).getLeft()+1),(snpInfo.get(parts2[1]).getRight()+1)));
                } else {
                    snpInfo.put(parts2[1], new Pair<>((snpInfo.get(parts2[1]).getLeft()+1),(snpInfo.get(parts2[1]).getRight())));
                }
            }

        }
        
        input.close();
        return new Pair<>(cpgInfo, snpInfo);
    }

    private static void processRealData(String inputFileReal, Pair<HashMap<String, Pair<Integer, Integer>>, HashMap<String, Pair<Integer, Integer>>> informationPerEntry, String outputFile) throws IOException {
        TextFile input = new TextFile(inputFileReal, TextFile.R);
        TextFile output = new TextFile(outputFile, TextFile.W);
        
        String row;
        while((row=input.readLine())!=null){
            String[] parts = row.split("\t");
            String[] parts2 = parts[0].split("-");
            
            output.writeln(row+"\t"+informationPerEntry.getLeft().get(parts2[0]).getLeft()+"\t"+informationPerEntry.getLeft().get(parts2[0]).getRight()+"\t"+informationPerEntry.getRight().get(parts2[1]).getLeft()+"\t"+informationPerEntry.getRight().get(parts2[1]).getRight());

        }
        output.close();
        input.close();
    }
    
    private static void processRealData2(Pair<HashMap<String, Pair<Integer, Integer>>, HashMap<String, Pair<Integer, Integer>>> informationPerRealEntry, Pair<HashMap<String, Pair<Integer, Integer>>, HashMap<String, Pair<Integer, Integer>>> informationPerPermEntry, String outputFile) throws IOException {
        TextFile outputSnp = new TextFile(outputFile+"snp.txt", TextFile.W);
        
        for(Entry<String, Pair<Integer, Integer>> snp : informationPerRealEntry.getRight().entrySet()){
            outputSnp.writeln(snp.getKey()+"\t"+snp.getValue().getLeft()+"\t"+snp.getValue().getRight()+"\t"+informationPerPermEntry.getRight().get(snp.getKey()).getLeft()+"\t"+informationPerPermEntry.getRight().get(snp.getKey()).getRight());
        }
        outputSnp.close();
        
        TextFile outputProbe = new TextFile(outputFile+"probe.txt", TextFile.W);
        for(Entry<String, Pair<Integer, Integer>> probe : informationPerRealEntry.getLeft().entrySet()){
            outputProbe.writeln(probe.getKey()+"\t"+probe.getValue().getLeft()+"\t"+probe.getValue().getRight()+"\t"+informationPerPermEntry.getLeft().get(probe.getKey()).getLeft()+"\t"+informationPerPermEntry.getLeft().get(probe.getKey()).getRight());
        }
        outputProbe.close();
        
    }
    
}
