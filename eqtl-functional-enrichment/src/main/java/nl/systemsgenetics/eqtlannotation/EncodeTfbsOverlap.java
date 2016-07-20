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
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.lang3.StringUtils;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;

/**
 *
 * @author MarcJan
 */
public class EncodeTfbsOverlap {

    public static void main(String[] args) {
        String inputFile = "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Trans_Pc22c_CisMeQTLc_meQTLs\\RegressedOut_CisEffects_New\\eQTLsFDR0.05-ProbeLevel_Filtered.txt";
        String outputFile = "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Trans_Pc22c_CisMeQTLc_meQTLs\\RegressedOut_CisEffects_New\\eQTLsFDR0.05-ProbeLevel_Filtered_NFKB_all_probeWindow.txt"
                + ".txt";
        String snpToCheck = "";
        String inputFolderTfbsData = "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\RP3_TFBS_NFKB1\\";
        int window = 25;
        HashMap<String, ArrayList<EncodeNarrowPeak>> peakData = null;

        ArrayList<EQTL> qtls = new ArrayList<>();
        try {
            QTLTextFile e = new QTLTextFile(inputFile, QTLTextFile.R);
            qtls.addAll(Arrays.asList(e.read()));
            e.close();

//            peakData = processOverlappingPeaks(readTfbsInformation(inputFolderTfbsData));
            peakData = readTfbsInformation(inputFolderTfbsData);

        } catch (IOException ex) {
            Logger.getLogger(EncodeTfbsOverlap.class.getName()).log(Level.SEVERE, null, ex);
        }

        try {
            TextFile outWriter = new TextFile(outputFile, TextFile.W);
            outWriter.writeln(QTLTextFile.header+"\tNFKB");
            for (EQTL e : qtls) {
//                System.out.println(e);
                if(e.getRsName().equals(snpToCheck) || snpToCheck.equals("")){
                    outWriter.writeln(e.toString() + "\t" + determineContact2(e, peakData, window));
                }
            }
            outWriter.close();
        } catch (IOException ex) {
            Logger.getLogger(EncodeTfbsOverlap.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    private static HashMap<String, ArrayList<EncodeNarrowPeak>> readTfbsInformation(String inputFolderTfbsData) throws IOException {
        HashMap<String, ArrayList<EncodeNarrowPeak>> data = new HashMap<>();
        File file = new File(inputFolderTfbsData);
        File[] files = file.listFiles();
        ArrayList<String> vecFiles = new ArrayList<>();
        for (File f : files) {
            System.out.println(f.getAbsolutePath());
            vecFiles.add(f.getAbsolutePath());
        }
        int totalSize = 0;
            
            for (String fileToRead : vecFiles) {
                TextFile reader = new TextFile(fileToRead, TextFile.R);
                
                String row;
                while((row=reader.readLine())!=null){
                    String[] parts = StringUtils.split(row, '\t');
                    
                    if(!data.containsKey(parts[0])){
                        data.put(parts[0], new ArrayList<EncodeNarrowPeak>());
                    }
                    data.get(parts[0]).add(new EncodeNarrowPeak(parts, fileToRead));
                    totalSize++;
                }
                
                reader.close();
        
        }
            System.out.println("Total contacts read in: "+totalSize);
        
        return data;
    }

    private static HashMap<String, ArrayList<EncodeNarrowPeak>> processOverlappingPeaks(HashMap<String, ArrayList<EncodeNarrowPeak>> peakData) {
        HashMap<String, ArrayList<EncodeNarrowPeak>> mergedPeakData = new HashMap<>();
        int totalMergedPeaks = 0;
        for(Entry<String, ArrayList<EncodeNarrowPeak>> chr : peakData.entrySet()){
            Collections.sort(chr.getValue());
            mergedPeakData.put(chr.getKey(), processPeaks(chr.getValue()));
            totalMergedPeaks+=mergedPeakData.get(chr.getKey()).size();
        }
        System.out.println("Total merged peaks: "+totalMergedPeaks);
        return(mergedPeakData);
        
    }

    private static ArrayList<EncodeNarrowPeak> processPeaks(ArrayList<EncodeNarrowPeak> value) {
        ArrayList<EncodeNarrowPeak> peaked = new ArrayList<>();
        
        for(EncodeNarrowPeak e : value){
            if(peaked.size()==0){
                peaked.add(e);
            } else {
              if(peaked.get((peaked.size()-1)).getChromEnd() <= e.getChromStart()){
                  peaked.set((peaked.size()-1), EncodeNarrowPeak.mergeTwoEntries(peaked.get((peaked.size()-1)), e));
              }  else {
                  peaked.add(e);
              }
            }
        }
        
        
        return peaked;
    }

    private static String determineContact(EQTL e, HashMap<String, ArrayList<EncodeNarrowPeak>> peakData) {
        String contact = "-\t-\t-";
        
        ArrayList<EncodeNarrowPeak> relevantData = peakData.get("chr"+e.getProbeChr());
        for(EncodeNarrowPeak peak : relevantData){
//            System.out.println(peak.getChromStart()+" "+peak.getChromEnd()+" "+e.getProbeChrPos());
            if(peak.getChromEnd()>= e.getProbeChrPos() && peak.getChromStart()<= e.getProbeChrPos()){
                contact = "overlapping\t"+peak.getName().split(",").length+"\t"+peak.getChrom()+"-"+peak.getChromStart()+"-"+peak.getChromEnd()+"\t"+peak.getName();
                break;
            } else if(peak.getChromEnd()> e.getProbeChrPos() && peak.getChromStart()> e.getProbeChrPos()){
                break;
            }
        }
        
        return contact;
    }
    
    private static String determineContact2(EQTL e, HashMap<String, ArrayList<EncodeNarrowPeak>> peakData, int window) {
        String contact = "-\t-\t-";
        
        ArrayList<String> name = new ArrayList<>();
        ArrayList<String> location = new ArrayList<>();
        
        ArrayList<EncodeNarrowPeak> relevantData = peakData.get("chr"+e.getProbeChr());
        for(EncodeNarrowPeak peak : relevantData){
//            System.out.println(peak.getChromStart()+" "+peak.getChromEnd()+" "+e.getProbeChrPos());
            if(peak.getChromEnd()>= (e.getProbeChrPos()-window) && peak.getChromStart()<= (e.getProbeChrPos()+window)){
                name.add(peak.getName());
                location.add(peak.getChrom()+"-"+peak.getChromStart()+"-"+peak.getChromEnd());
            }
        }
        if(!name.isEmpty()){
            StringBuilder allNames = new StringBuilder(name.get(0));
            StringBuilder allLocations = new StringBuilder(location.get(0));
            for(int i=1; i<name.size(); i++){
                allNames.append(",").append(name.get(i));
                allLocations.append(",").append(location.get(i));
            }
            contact = "overlapping\t"+name.size()+"\t"+allLocations.toString()+"\t"+allNames.toString();
        }
        
        
        return contact;
    }
}
