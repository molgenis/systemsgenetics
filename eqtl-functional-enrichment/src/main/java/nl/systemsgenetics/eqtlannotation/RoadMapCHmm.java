/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlannotation;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author MarcJan
 */
public class RoadMapCHmm {
    static final Pattern TAB_PATTERN = Pattern.compile("\\t");
    
    public static void main(String[] args) {
        String inputFile = "D:\\T\\eQTL\\eQTLs_with_significant_modules.txt";
        String outputFile = "D:\\T\\eQTL\\eQTLs_with_significant_modules_Annotated_S25_Blood";
        String inputFolderTfbsData = "D:\\T\\eQTL\\25S\\blood.mnemonics.bedFiles\\";
        
        HashMap<String, LinkedHashMap<String, ArrayList<ChmmState>>>  peakData = null;

        ArrayList<SNP> qtlSnps = new ArrayList<>();
        try {
            
            qtlSnps.addAll(readDataFromQTLFile(inputFile));

            peakData = readCHmmInformation(inputFolderTfbsData);

        } catch (IOException ex) {
            Logger.getLogger(EncodeTfbsOverlap.class.getName()).log(Level.SEVERE, null, ex);
        }

//        writeContacts(qtlSnps, peakData, outputFile);
        writeContacts2(qtlSnps, peakData, outputFile);

    }
    
    private static ArrayList<SNP> readDataFromQTLFile(String inputFile) throws IOException {
        ArrayList<SNP> qtlSnp = new ArrayList<>();
        TextFile f = new TextFile(inputFile, TextFile.R);
        
        String s=f.readLine();
        
        while((s=f.readLine())!=null){
//            System.out.println(s);
            qtlSnp.add(new SNP(TAB_PATTERN.split(s)));
        }
        
        
        f.close();
        return qtlSnp;
    }

    private static HashMap<String, LinkedHashMap<String, ArrayList<ChmmState>>> readCHmmInformation(String inputFolderTfbsData) throws IOException {
        //  chr, chmmState-tissue
        HashMap<String, LinkedHashMap<String, ArrayList<ChmmState>>> data = new HashMap<>();
        File file = new File(inputFolderTfbsData);
        File[] files = file.listFiles();
        
        ArrayList<String> vecFiles = new ArrayList<>();
        for (File f : files) {
            System.out.println(f.getName());
            vecFiles.add(f.getName());
        }
        
        int totalSize = 0;
//            
        for (String fileToRead : vecFiles) {
            TextFile reader = new TextFile(inputFolderTfbsData+File.separator+fileToRead, TextFile.R);

            String row;
            while((row=reader.readLine())!=null){
//                System.out.println(row);
                String[] parts = TAB_PATTERN.split(row);
                
                if(!data.containsKey(parts[0])){
                    LinkedHashMap<String, ArrayList<ChmmState>> perChr = new LinkedHashMap<>();
                    data.put(parts[0], perChr);
                    
                }
                if(!data.get(parts[0]).containsKey(parts[3]+"-"+fileToRead)){
                    ArrayList<ChmmState> perTissueState = new ArrayList<>();
                    data.get(parts[0]).put((parts[3]+"-"+fileToRead), perTissueState);
                }
                
                
                data.get(parts[0]).get(parts[3]+"-"+fileToRead).add(new ChmmState(parts));
                totalSize++;
            }
            reader.close();
        
        }
            System.out.println("Total contacts read in: "+totalSize);
        
        return data;
    }

    private static void determineContact(long pos, ArrayList<ChmmState> peakData, StringBuilder row, String key) {
        
        for(ChmmState peak : peakData){
            if(peak.getStop()>= pos && peak.getStart()<= pos){
                row.append('\t').append(key).append(':').append("+");
                break;
            } else if(peak.getStop()> pos && peak.getStart()> pos){
                break;
            }
        }
    }

    private static int determineContact(long pos, ArrayList<ChmmState> peakData) {
        
        for(ChmmState peak : peakData){
            if(peak.getStop()>= pos && peak.getStart()<= pos){
                return 1;
            } else if(peak.getStop()> pos && peak.getStart()> pos){
                break;
            }
        }
        return 0;
    }
    
    private static void writeContacts(ArrayList<SNP> qtlSnps, HashMap<String, LinkedHashMap<String, ArrayList<ChmmState>>> peakData, String outputFile) {
        try {
            TextFile outWriter = new TextFile(outputFile, TextFile.W);
            outWriter.writeln();
            for (SNP s : qtlSnps) {
                StringBuilder row = new StringBuilder();
                row.append(s.getPrimaryId()).append('\t');
                
                int i = 0;
                row.append(s.getDirection()[i]).append(s.getModuleAssociation()[i]);
                i++;
                for(; i<s.getModuleAssociation().length; ++i){
                    row.append(',').append(s.getDirection()[i]).append(s.getModuleAssociation()[i]);
                }
//                for(int module : s.getModuleAssociation()){
//                    //                System.out.println(e);
                for(Entry<String, ArrayList<ChmmState>> data : peakData.get("chr"+s.getChr()).entrySet()){
                    determineContact(s.getPosition(), data.getValue(), row,data.getKey());
                }
                row.append('\n');
//                }
//                System.out.print(row.toString());
                outWriter.write(row.toString());
            }
            outWriter.close();
        } catch (IOException ex) {
            Logger.getLogger(EncodeTfbsOverlap.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    private static void writeContacts2(ArrayList<SNP> qtlSnps, HashMap<String, LinkedHashMap<String, ArrayList<ChmmState>>> peakData, String outputFile) {
        HashMap<String,HashMap<String, HashMap<String, Integer>>> perChmmInformation = new HashMap<>();
        //chmm state // tissue //module // count.
        HashMap<String, HashMap<String, Integer>> BackGroundInformation = new HashMap<>();
        
        for (SNP s : qtlSnps) {
            for(Entry<String, ArrayList<ChmmState>> data : peakData.get("chr"+s.getChr()).entrySet()){
//                System.out.println(data.getKey());
                int overlap = determineContact(s.getPosition(), data.getValue());
                
                String[] parts = data.getKey().split("-");
                String chmm = parts[0];
                chmm = chmm.replace("/", "_");
                String tissue = parts[1].split("_")[0];
                
                if(!BackGroundInformation.containsKey(tissue)){
                    BackGroundInformation.put(tissue, new HashMap<String, Integer>());
                }
                
                if(!perChmmInformation.containsKey(chmm)){
                    perChmmInformation.put(chmm, new HashMap<String, HashMap<String, Integer>>());
                }
                for(int i = 0; i<s.getModuleAssociation().length; ++i){
                    String Module = String.valueOf(s.getModuleAssociation()[i]);
                    if(!BackGroundInformation.get(tissue).containsKey(Module)){
                        BackGroundInformation.get(tissue).put(Module, 0);
                    }
                    BackGroundInformation.get(tissue).put(Module,(BackGroundInformation.get(tissue).get(Module)+overlap));
                    
                    if(!perChmmInformation.get(chmm).containsKey(tissue)){
                        perChmmInformation.get(chmm).put(tissue, new HashMap<String, Integer>());
                    }
                    if(!perChmmInformation.get(chmm).get(tissue).containsKey(Module)){
                        perChmmInformation.get(chmm).get(tissue).put(Module, 0);
                    }
                    perChmmInformation.get(chmm).get(tissue).put(Module, (perChmmInformation.get(chmm).get(tissue).get(Module)+overlap));
                }
            }
        }

        String[] modules = {"0","1","2","3","4","5","6","7","8","9","10"};
        try{
            LinkedHashSet<String> keys = new LinkedHashSet<>();
            for(Entry<String, HashMap<String,Integer>> e :BackGroundInformation.entrySet()){
                keys.addAll(e.getValue().keySet());
            }
            
            TextFile outWriter = new TextFile(outputFile+"_Background.txt", TextFile.W);
            StringBuilder s = new StringBuilder();
            s.append("Tissue");
            
            for(String k : modules){
                s.append("\t").append(k);
            }
            outWriter.writeln(s.toString());
            
            for(Entry<String, HashMap<String,Integer>> e :BackGroundInformation.entrySet()){
                s = new StringBuilder();
                s.append(e.getKey());
                for(String k : modules){
                    if(e.getValue().containsKey(k)){
                        s.append('\t').append(e.getValue().get(k));
                    } else {
                        s.append('\t').append(0);
                    }
                }
                outWriter.writeln(s.toString());
            }
            outWriter.close();
        } catch (IOException ex) {
            Logger.getLogger(EncodeTfbsOverlap.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        try{
            //chmm state //module // tissue // count.
            for(Entry<String,HashMap<String, HashMap<String, Integer>>> e: perChmmInformation.entrySet()){
                TextFile outWriter = new TextFile(outputFile+"_"+e.getKey()+".txt", TextFile.W);
                
                StringBuilder s = new StringBuilder();
                s.append("Tissue");

                for(String k : modules){
                    s.append("\t").append(k);
                }
                outWriter.writeln(s.toString());

                for(Entry<String, HashMap<String,Integer>> f :e.getValue().entrySet()){
                    s = new StringBuilder();
                    s.append(f.getKey());
                    for(String k : modules){
                        if(f.getValue().containsKey(k)){
                            s.append('\t').append(f.getValue().get(k));
                        } else {
                            s.append('\t').append(0);
                        }
                    }
                    outWriter.writeln(s.toString());
                }
                outWriter.close();
            }
        } catch (IOException ex) {
            Logger.getLogger(EncodeTfbsOverlap.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

}
