/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlannotation;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.lang3.StringUtils;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;
import umcg.genetica.io.trityper.util.ChrAnnotation;

/**
 *
 * @author MarcJan
 */
public class EncodeMultipleTfbsOverlap {

    public static void main(String[] args) {
//        String inputFile = "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Trans_Pc22c_CisMeQTLc_meQTLs\\RegressedOut_CisEffects_New\\Permutations\\PermutedEQTLsPermutationRoundX.head.extended.txt";
//        String outputFile = "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Trans_Pc22c_CisMeQTLc_meQTLs\\RegressedOut_CisEffects_New\\Permutations\\PermutedEQTLsPermutationRoundX.head.extended.txtall_probeWindow_Test.txt";
        String inputFile = "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Trans_Pc22c_CisMeQTLc_meQTLs\\RegressedOut_CisEffects_New\\eQTLsFDR0.05-ProbeLevel_Filtered.txt";
        String outputFile = "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Trans_Pc22c_CisMeQTLc_meQTLs\\RegressedOut_CisEffects_New\\TFBS\\Test.txt";
        
//        String inputFile = "E:\\LL+RS+CODAM+LLS_eqtls_genes_allLevels_18122014_LD_99_TF_input.txt";
//        String outputFile = "E:\\LL+RS+CODAM+LLS_eqtls_genes_allLevels_18122014_LD_99_TF_input_all.txt";
        
        
        String snpToCheck = "";
//        String inputFolderTfbsData = "E:\\OnlineFolders\\BitSync\\UMCG\\Data\\ENCODE_TFBS\\";
        String inputFolderTfbsData = "D:\\UMCG\\Data\\ENCODE_TFBS\\";
        int window = 1;
        
        boolean perm = false;
        
        LinkedHashMap<String,HashMap<String, ArrayList<EncodeNarrowPeak>>> peakData = null;

        ArrayList<EQTL> qtls = new ArrayList<>();
        try {
            if(perm){
                TextFile t = new TextFile(inputFile, TextFile.R);
                String row = t.readLine();
                while((row=t.readLine())!=null){
                    String[] parts = row.split("\t");
                    EQTL e = new EQTL();
//                    System.out.println(row);
                    e.setRsName(parts[1]);
                    e.setRsChr(ChrAnnotation.parseChr(parts[2]));
                    e.setRsChrPos(Integer.parseInt(parts[3]));
                    e.setProbe(parts[4]);
                    e.setProbeChr(ChrAnnotation.parseChr(parts[5]));
                    e.setProbeChrPos(Integer.parseInt(parts[6]));
                    
                    qtls.add(e);
                }
                
                t.close();
            } else {
                QTLTextFile e = new QTLTextFile(inputFile, QTLTextFile.R);
                qtls.addAll(Arrays.asList(e.read()));
                e.close();
            }

//            peakData = processOverlappingPeaks(readTfbsInformation(inputFolderTfbsData));
            peakData = readMultipleTfbsInformation(inputFolderTfbsData);

        } catch (IOException ex) {
            Logger.getLogger(EncodeMultipleTfbsOverlap.class.getName()).log(Level.SEVERE, null, ex);
        }

        try {
            TextFile outWriter = new TextFile(outputFile, TextFile.W);
            
            StringBuilder extraHeaderInfo = new StringBuilder();
            for(String keys : peakData.keySet()){
                extraHeaderInfo.append('\t').append(keys).append('\t').append("totalOverlaps").append('\t').append("bindingLocations").append('\t').append("fileNames");
            }
            
            outWriter.writeln(QTLTextFile.header+extraHeaderInfo.toString());
            for (EQTL e : qtls) {
//                System.out.println(e);
                if(e.getRsName().equals(snpToCheck) || snpToCheck.equals("")){
                    outWriter.writeln(e.toString() + determineContacts(e, peakData, window));
                }
            }
            outWriter.close();
        } catch (IOException ex) {
            Logger.getLogger(EncodeMultipleTfbsOverlap.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    private static LinkedHashMap<String,HashMap<String, ArrayList<EncodeNarrowPeak>>> readMultipleTfbsInformation(String inputFolderTfbsData) throws IOException {
        LinkedHashMap<String,HashMap<String, ArrayList<EncodeNarrowPeak>>> data = new LinkedHashMap<>();
        File file = new File(inputFolderTfbsData);
        File[] files = file.listFiles();
        ArrayList<String> vecFiles = new ArrayList<>();
        for (File f : files) {
//            System.out.println(f.getAbsolutePath());
            vecFiles.add(f.getAbsolutePath());
        }
        
        for (String fileToRead : vecFiles) {
            TextFile reader = new TextFile(fileToRead, TextFile.R);
            
            String[] storingInformation = fileToRead.split("_");
//            String cellLine = storingInformation[1].replace("TFBS\\","");
            String transcriptionFactor = storingInformation[2].replace(".narrowPeak","");
            if(storingInformation.length>4){
                for(int i=3;i<(storingInformation.length-1);++i){
                    transcriptionFactor = transcriptionFactor +"_"+storingInformation[i].replace(".narrowPeak","");
                }
            }
            
            String row;
            while((row=reader.readLine())!=null){

                String[] parts = StringUtils.split(row, '\t');
                if(!data.containsKey(transcriptionFactor)){
                    data.put(transcriptionFactor, new HashMap<String, ArrayList<EncodeNarrowPeak>>());
                }
                if(!data.get(transcriptionFactor).containsKey(parts[0])){
                    data.get(transcriptionFactor).put(parts[0], new ArrayList<EncodeNarrowPeak>());
                }
                data.get(transcriptionFactor).get(parts[0]).add(new EncodeNarrowPeak(parts, fileToRead));
            }

            reader.close();
        
        }
        ArrayList<String> cleanList = new ArrayList<>();
        for(Entry<String,HashMap<String, ArrayList<EncodeNarrowPeak>>> tfInformation : data.entrySet()){
            System.out.println("Transcription factor: "+tfInformation.getKey());
            int counter = 0;
            for(Entry<String, ArrayList<EncodeNarrowPeak>> tfEntry : tfInformation.getValue().entrySet()){
                Collections.sort(tfEntry.getValue());
                counter+=tfEntry.getValue().size();
            }
            System.out.println("\tcontacts: "+counter);
            
            //remove all with less than 750 contacts
//            if(counter<750){
//                cleanList.add(tfInformation.getKey());
//            }
        }
        
        for(String k : cleanList){
            data.remove(k);
        }
        
        return data;
    }

    private static String determineContacts(EQTL e, HashMap<String,HashMap<String, ArrayList<EncodeNarrowPeak>>> peakData, int window) {
        StringBuilder returnableContacts = new StringBuilder();
        
        for(Entry<String,HashMap<String, ArrayList<EncodeNarrowPeak>>> tfData : peakData.entrySet()){
            returnableContacts.append('\t');
            ArrayList<String> name = new ArrayList<>();
            ArrayList<String> location = new ArrayList<>();
            if(tfData.getValue().containsKey("chr"+e.getProbeChr())){
                ArrayList<EncodeNarrowPeak> relevantData = tfData.getValue().get("chr"+e.getProbeChr());
                if(!relevantData.isEmpty()){
                    for(EncodeNarrowPeak peak : relevantData){
            //            System.out.println(peak.getChromStart()+" "+peak.getChromEnd()+" "+e.getProbeChrPos());
                        if((peak.getChromEnd()+window)>= e.getProbeChrPos() && (peak.getChromStart()-window)<= (e.getProbeChrPos())){
                            name.add(peak.getName());
                            location.add(peak.getChrom()+"-"+peak.getChromStart()+"-"+peak.getChromEnd());
                        } else if((peak.getChromEnd()+window)> e.getProbeChrPos() && (peak.getChromStart()-window)> e.getProbeChrPos()){
                            break;
                        }
                    }
                }
                if(!name.isEmpty()){
                    StringBuilder allNames = new StringBuilder(name.get(0));
                    StringBuilder allLocations = new StringBuilder(location.get(0));
                    for(int i=1; i<name.size(); i++){
                        allNames.append(",").append(name.get(i));
                        allLocations.append(",").append(location.get(i));
                    }
                    returnableContacts.append("overlapping\t").append(name.size()).append('\t').append(allLocations.toString()).append('\t').append(allNames.toString());
                } else {
                    returnableContacts.append("\t\t\t");
                }
            } else {
                returnableContacts.append("\t\t\t");
            }
        }
        return returnableContacts.toString();
    }
}
