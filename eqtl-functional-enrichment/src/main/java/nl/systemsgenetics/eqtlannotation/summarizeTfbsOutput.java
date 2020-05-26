/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlannotation;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.FisherExactTest;

/**
 *
 * @author MarcJan
 */
public class summarizeTfbsOutput {

    public static void main(String[] args) {
        
        String inputFile = "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Trans_Pc22c_CisMeQTLc_meQTLs\\RegressedOut_CisEffects_New\\TFBS\\rs11190133\\eQTLsFDR0.05-ProbeLevel_Filtered_rs11190133_Interest.txt";
        String outputFile = "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Trans_Pc22c_CisMeQTLc_meQTLs\\RegressedOut_CisEffects_New\\TFBS\\rs11190133\\AllFactor_Summary.txt";
        
        try {
            summarize(inputFile, outputFile);
        } catch (IOException ex) {
            Logger.getLogger(summarizeTfbsOutput.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    private static void summarize(String inputFile, String outputFile) throws IOException {
        TextFile reader = new TextFile(inputFile, TextFile.R);
        String info = reader.readLine();
        String[] headerParts = info.split("\t");
        
        HashSet<String> allCpGs = new HashSet<>();
        HashMap<String, HashSet<String>> qtls = new HashMap<>();
        
        HashMap<String,HashSet<String>> cpgStats = new HashMap<>();

        while ((info = reader.readLine()) != null) {
            String[] parts = info.split("\t");
            
            String snpName = parts[1];
            String cpGName = parts[4];
            allCpGs.add(cpGName);
            
            
            if(!qtls.containsKey(snpName)){
                qtls.put(snpName, new HashSet<String>());
            } 
            qtls.get(snpName).add(cpGName);
            
            if(!cpgStats.containsKey(cpGName)){
                cpgStats.put(cpGName, new HashSet<String>());
                for (int i = 22; (i < headerParts.length && i < parts.length); i += 4) {
                    
                    if(parts[i].equals("overlapping")){
                        cpgStats.get(cpGName).add(headerParts[i]);
                    }
                }
            }
            
            
            String dir = "-";
            double zScore = Double.parseDouble(parts[10]);
            if(zScore>0){
                dir = "+";
            }
            snpName = parts[1]+dir;
            
            if(!qtls.containsKey(snpName)){
                qtls.put(snpName, new HashSet<String>());
            } 
            qtls.get(snpName).add(cpGName);
            
            if(!cpgStats.containsKey(cpGName)){
                cpgStats.put(cpGName, new HashSet<String>());
                for (int i = 22; (i < headerParts.length && i < parts.length); i += 4) {
                    
                    if(parts[i].equals("overlapping")){
                        cpgStats.get(cpGName).add(headerParts[i]);
                    }
                }
            }
            
            
        }
        reader.close();
        
        writeSummary(outputFile, qtls, allCpGs, cpgStats);
        
    }

    private static void writeSummary(String outputFile, HashMap<String, HashSet<String>> qtls, HashSet<String> allCpGs, HashMap<String, HashSet<String>> cpgStats) throws IOException {
        
        TextFile writer = new TextFile(outputFile, TextFile.W);
        FisherExactTest f = new FisherExactTest();
        
        for(Entry<String, HashSet<String>> qtlSnp : qtls.entrySet()){
            ArrayList<TfbsOverlapResult> resultBuffer = new ArrayList<>();
            
            HashSet<String> traits = new HashSet<>();
            for(String CpG : qtlSnp.getValue()){
                traits.addAll(cpgStats.get(CpG));
            }
            
            //Loop over all relevant traits
            for(String trait : traits){
                int nrOfTimesTfForSnp = 0;
                int nrOfTimesSnpNoneTf = 0;
                int nrOfTimesTfFOtherSnps = 0;
                int currentOther = 0;
                //Loop over all CpGs.
                for(String cpg : allCpGs){
                    if(qtlSnp.getValue().contains(cpg)){
                        if(cpgStats.get(cpg).contains(trait)){
                            nrOfTimesTfForSnp++;
                        } else {
                            nrOfTimesSnpNoneTf++;
                        }
                    } else {
                        if(cpgStats.get(cpg).contains(trait)){
                            nrOfTimesTfFOtherSnps++;
                        } else {
                            currentOther++;
                        }
                    }
                }
                
                double fisher = f.getFisherPValue(nrOfTimesTfForSnp,nrOfTimesSnpNoneTf,nrOfTimesTfFOtherSnps,currentOther);
//                System.out.println(new TfbsOverlapResult(trait,fisher,nrOfTimesTfForSnp, nrOfTimesSnpNoneTf, nrOfTimesTfFOtherSnps, currentOther).toString());
                if((nrOfTimesTfForSnp+nrOfTimesSnpNoneTf)>10 &&fisher<=0.05){
                    resultBuffer.add(new TfbsOverlapResult(trait,fisher,nrOfTimesTfForSnp, nrOfTimesSnpNoneTf, nrOfTimesTfFOtherSnps, currentOther));
                }
            }

            if(resultBuffer.size()>0){
                Collections.sort(resultBuffer);
                writer.writeln(qtlSnp.getValue().size()+"\t"+qtlSnp.getKey());
                System.out.println(qtlSnp.getKey());
                for(TfbsOverlapResult r : resultBuffer){
                    writer.writeln("\t\t\t"+r.toString());
                    System.out.println("\t\t\t"+r.toString());
                }
            }
        }
        writer.close();
    }

}
