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
        String inputFile = "F:\\eQTL\\eQTLs_with_significant_modules.txt";
        String outputFile = "F:\\eQTL\\eQTLs_with_significant_modules_Annotated_S15_Blood.txt";
        String inputFolderTfbsData = "F:\\eQTL\\15S\\blood.mnemonics.bedFiles\\";
        
        HashMap<String, LinkedHashMap<String, ArrayList<ChmmState>>>  peakData = null;

        ArrayList<SNP> qtlSnps = new ArrayList<>();
        try {
            
            qtlSnps.addAll(readDataFromQTLFile(inputFile));

            peakData = readCHmmInformation(inputFolderTfbsData);

        } catch (IOException ex) {
            Logger.getLogger(EncodeTfbsOverlap.class.getName()).log(Level.SEVERE, null, ex);
        }

        writeContacts(qtlSnps, peakData, outputFile);

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


}
