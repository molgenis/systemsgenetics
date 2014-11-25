/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.util;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author MarcJan
 */
public class QTLLookup {
    
    public static void main(String[] args) throws IOException {
        lookUpEffects(
                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\QTLCorrected\\RP3_2MB_TSS_extendedCis_eQTMs\\eQTLsFDR-SNPLevel.txt.gz",
                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\DEPICTAnalysisrs3774959-250kbAroundCpGAnalysis_CpG-ProbeList.txt",
                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\DEPICT_eQTM_Strengts_rs3774959_cisCorrected.txt");

    }

    static void lookUpEffects(String in, String in2, String out) throws IOException {
        TextFile inCombi = new TextFile(in2, TextFile.R);
        HashSet<String> combiList = new HashSet<String>(inCombi.readAsArrayList());
        
        HashSet<String> identified = new HashSet<String>();
        
        TextFile e = new TextFile(in, TextFile.R);

        TextFile outWriter = new TextFile(out, TextFile.W);
        outWriter.writeln("Combination\tP-value\tZscore\tFDR\tStrongesteQTMOnProbe\tP-value\tZscore\tFDR\tSameEffect");
        
        HashMap<String, String> strongestEffectPerProbe = new HashMap<String, String>();
        String str;
        while((str = e.readLine())!=null){
            String[] parts = str.split(("\t"));
            if(!strongestEffectPerProbe.containsKey(parts[1])){
                strongestEffectPerProbe.put(parts[1], (parts[4]+"\t"+parts[0]+"\t"+parts[10]+"\t"+parts[(parts.length-1)]));
            }
            
            if(combiList.contains(parts[1]+"-"+parts[4])){
                String partToWrite = parts[4]+"\t"+parts[0]+"\t"+parts[10]+"\t"+parts[(parts.length-1)];
                if(partToWrite.equals(strongestEffectPerProbe.get(parts[1]))){
                    outWriter.writeln(parts[1]+"-"+partToWrite+"\t"+strongestEffectPerProbe.get(parts[1])+"\tIdentical");
                } else {
                    outWriter.writeln(parts[1]+"-"+partToWrite+"\t"+strongestEffectPerProbe.get(parts[1])+"\tnonIdentical");
                }
                identified.add(parts[1]+"-"+parts[4]);
            }
        }
        
        for(String combination : combiList){
            if(!identified.contains(combination)&&strongestEffectPerProbe.containsKey(combination.split("-")[0])){
                outWriter.writeln(combination+"\tNA\tNA\tNA\t"+strongestEffectPerProbe.get(combination.split("-")[0])+"\t-");
            } else if (!identified.contains(combination)) {
                outWriter.writeln(combination+"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t-");
            }
        }
        outWriter.close();
    }
    
}
