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
//        lookUpEffects(
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\QTLCorrected\\RP3_2MB_TSS_extendedCis_eQTMs\\eQTLsFDR-SNPLevel.txt.gz",
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\DEPICTAnalysisrs3774959-250kbAroundCpGAnalysis_CpG-ProbeList.txt",
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\DEPICT_eQTM_Strengts_rs3774959_cisCorrected.txt");
        
//        lookUpEffects2(
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\QTLCorrected\\RP3_2MB_TSS_extendedCis_eQTMs\\eQTLsFDR-SNPLevel.txt.gz",
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\Great_GeneAssociations_20141127-public-2.0.2-aaQ50i-hg19-all-region_Ensembl.txt",
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\Great_eQTM_Strengts_cisCorrected2mb.txt");
        
        lookUpEffects3(
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Cis_Pc22c_meQTLs\\eQTLProbesFDR0.05-1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-17-ProbeLevel_noLD.txt",
                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Cis_Pc22c_meQTLs\\eQTLProbesFDR0.05-1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-ProbeLevel.txt",
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Trans_Pc22c_CisMeQTLc_meQTLs\\RegressedOut_CisEffects_New\\eQTLsFDR0.05-ProbeLevel.txt",
                "D:\\UMCG\\Requests\\EWAS\\Files_for_MarcJan\\Genetic_profile\\BMI_2014.txt",
//                null,
                "D:\\UMCG\\Requests\\EWAS\\Files_for_MarcJan\\OutputDirSonia\\probeList.txt",
                "D:\\UMCG\\Requests\\EWAS\\Files_for_MarcJan\\OutputDirSonia\\LookedUpEffects_Cis2.txt");
    }

    static void lookUpEffects(String in, String in2, String out) throws IOException {
        TextFile inCombi = new TextFile(in2, TextFile.R);
        HashSet<String> combiList = new HashSet<String>(inCombi.readAsArrayList());
        
        HashSet<String> identified = new HashSet<String>();
        
        TextFile e = new TextFile(in, TextFile.R);

        TextFile outWriter = new TextFile(out, TextFile.W);
        outWriter.writeln("CpG\tProbe\tCombination\tP-value\tZscore\tFDR\tStrongesteQTMOnProbe\tP-value\tZscore\tFDR\tSameEffect");
        
        HashMap<String, String> strongestEffectPerProbe = new HashMap<String, String>();
        HashSet<String> allCpGs = new HashSet<String>();
        HashSet<String> identicalCpGs = new HashSet<String>();
                
        String str;
        while((str = e.readLine())!=null){
            String[] parts = str.split(("\t"));
            allCpGs.add(parts[1]);
            if(!strongestEffectPerProbe.containsKey(parts[1])){
                strongestEffectPerProbe.put(parts[1], (parts[4]+"\t"+parts[0]+"\t"+parts[10]+"\t"+parts[(parts.length-1)]));
            }
            
            if(combiList.contains(parts[1]+"-"+parts[4])){
                String partToWrite = parts[4]+"\t"+parts[0]+"\t"+parts[10]+"\t"+parts[(parts.length-1)];
                if(partToWrite.equals(strongestEffectPerProbe.get(parts[1]))){
                    outWriter.writeln(parts[1]+"\t"+parts[4]+"\t"+parts[1]+"-"+partToWrite+"\t"+strongestEffectPerProbe.get(parts[1])+"\tIdentical");
                    identicalCpGs.add(parts[1]);
                } else {
                    outWriter.writeln(parts[1]+"\t"+parts[4]+"\t"+parts[1]+"-"+partToWrite+"\t"+strongestEffectPerProbe.get(parts[1])+"\tnonIdentical");
                }
                identified.add(parts[1]+"-"+parts[4]);
            }
        }
        
        for(String combination : combiList){
            String[] parts = combination.split("-");
            if(!identified.contains(combination)&&strongestEffectPerProbe.containsKey(combination.split("-")[0])){
                outWriter.writeln(parts[0]+"\t"+parts[1]+"\t"+combination+"\tNA\tNA\tNA\t"+strongestEffectPerProbe.get(combination.split("-")[0])+"\t-");
            } else if (!identified.contains(combination)) {
                outWriter.writeln(parts[0]+"\t"+parts[1]+"\t"+"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t-");
            }
        }
        
        System.out.println((identicalCpGs.size()/(double)allCpGs.size())*100);
        
        outWriter.close();
    }
    
    static void lookUpEffects2(String in, String in2, String out) throws IOException {
        TextFile inCombi = new TextFile(in2, TextFile.R);
        HashSet<String> tmpCombiList = new HashSet<String>(inCombi.readAsArrayList());
        HashSet<String> combiList = new HashSet<String>();
        for(String s : tmpCombiList){
            if(s.contains("\t") && s.contains("ENSG")){
                String[] parts = s.split("\t");
                if(parts[1].contains(", ")){
                    
                    String subparts[] = parts[1].split(", ");
                    for(String t : subparts){
                        combiList.add(parts[0]+"-"+t);
//                        System.out.println(parts[0]+"-"+t);
                    }
                } else {
                    combiList.add(parts[0]+"-"+parts[1]);
//                    System.out.println(parts[0]+"-"+parts[1]);
                }
            }
        }
        
        HashSet<String> identified = new HashSet<String>();
        HashSet<String> allCpGs = new HashSet<String>();
        HashSet<String> identicalCpGs = new HashSet<String>();
        TextFile e = new TextFile(in, TextFile.R);

        TextFile outWriter = new TextFile(out, TextFile.W);
        outWriter.writeln("CpG\tProbe\tCombination\tP-value\tZscore\tFDR\tStrongesteQTMOnProbe\tP-value\tZscore\tFDR\tSameEffect");
        
        HashMap<String, String> strongestEffectPerProbe = new HashMap<String, String>();
        HashMap<String, Double> strongestZscorePerProbe = new HashMap<String, Double>();
        String str = e.readLine();
        while((str = e.readLine())!=null){
            String[] parts = str.split(("\t"));
            if(!strongestEffectPerProbe.containsKey(parts[1])){
                strongestEffectPerProbe.put(parts[1], (parts[4]+"\t"+parts[0]+"\t"+parts[10]+"\t"+parts[(parts.length-1)]));
                strongestZscorePerProbe.put(parts[1], Double.parseDouble(parts[10]));
            }
            
            if(combiList.contains(parts[1]+"-"+parts[4])){
                allCpGs.add(parts[1]);
                String partToWrite = parts[4]+"\t"+parts[0]+"\t"+parts[10]+"\t"+parts[(parts.length-1)];
                if(partToWrite.equals(strongestEffectPerProbe.get(parts[1]))){
                    outWriter.writeln(parts[1]+"\t"+parts[4]+"\t"+parts[1]+"-"+partToWrite+"\t"+strongestEffectPerProbe.get(parts[1])+"\tIdentical\t0");
                    identicalCpGs.add(parts[1]);
                } else {
                    outWriter.writeln(parts[1]+"\t"+parts[4]+"\t"+parts[1]+"-"+partToWrite+"\t"+strongestEffectPerProbe.get(parts[1])+"\tnonIdentical\t"+Math.abs((Math.abs(Double.parseDouble(parts[10])) - Math.abs(strongestZscorePerProbe.get(parts[1])))));
                }
                identified.add(parts[1]+"-"+parts[4]);
            }
        }
        
        for(String combination : combiList){
            String[] parts = combination.split("-");
            allCpGs.add(parts[0]);
            if(!identified.contains(combination)&&strongestEffectPerProbe.containsKey(combination.split("-")[0])){
                outWriter.writeln(parts[0]+"\t"+parts[1]+"\t"+combination+"\tNA\tNA\tNA\t"+strongestEffectPerProbe.get(combination.split("-")[0])+"\t-");
            } else if (!identified.contains(combination)) {
                outWriter.writeln(parts[0]+"\t"+parts[1]+"\t"+combination+"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t-");
            }
        }
        outWriter.close();
        System.out.println(identicalCpGs.size());
        System.out.println((double)allCpGs.size());
    }
    
    static void lookUpEffects3(String inQTLs, String inLookupSnps, String inLookupCpGs, String out) throws IOException {
        if(inLookupSnps != null &&  inLookupCpGs != null){
            TextFile inSnps = new TextFile(inLookupSnps, TextFile.R);
            TextFile inProbes = new TextFile(inLookupCpGs, TextFile.R);

            HashSet<String> tmpSnpsList = new HashSet<String>(inSnps.readAsArrayList());
            HashSet<String> tmpProbesList = new HashSet<String>(inProbes.readAsArrayList());
            
            HashSet<String> newTmpSnpsList = new HashSet<String>();
            for(String s : tmpSnpsList){
                String newS = s;
                if(s.contains("\t")){
                    newS = s.split("\t")[0];
                }
                newTmpSnpsList.add(newS);
            }
            tmpSnpsList = null;
            
            
            TextFile e = new TextFile(inQTLs, TextFile.R);
            
            TextFile outWriter = new TextFile(out, TextFile.W);
            outWriter.writeln("SNP\tProbe\tCombination\tP-value\tZscore\tFDR");

            String str = e.readLine();
            while((str = e.readLine())!=null){
                String[] parts = str.split(("\t"));

                if(newTmpSnpsList.contains(parts[1]) && tmpProbesList.contains(parts[4])){
                    outWriter.writeln(parts[1]+"\t"+parts[4]+"\t"+parts[1]+"-"+parts[4]+"\t"+parts[0]+"\t"+parts[10]+"\t"+parts[(parts.length-1)]);
                }
            }

            outWriter.close();
        } else if (inLookupSnps != null) {
            TextFile inSnps = new TextFile(inLookupSnps, TextFile.R);

            HashSet<String> tmpSnpsList = new HashSet<String>(inSnps.readAsArrayList());
            HashSet<String> newTmpSnpsList = new HashSet<String>();
            for(String s : tmpSnpsList){
                String newS = s;
                if(s.contains("\t")){
                    newS = s.split("\t")[0];
                }
                newTmpSnpsList.add(newS);
            }
            tmpSnpsList = null;
            
            TextFile e = new TextFile(inQTLs, TextFile.R);

            TextFile outWriter = new TextFile(out, TextFile.W);
            outWriter.writeln("SNP\tCombination\tP-value\tZscore\tFDR");

            String str = e.readLine();
            while((str = e.readLine())!=null){
                String[] parts = str.split(("\t"));

                if(newTmpSnpsList.contains(parts[1])){
                    outWriter.writeln(parts[1]+"\t"+parts[1]+"-"+parts[4]+"\t"+parts[0]+"\t"+parts[10]+"\t"+parts[(parts.length-1)]);
                }
            }

            outWriter.close();
        } else if(inLookupCpGs != null){
            TextFile inProbes = new TextFile(inLookupCpGs, TextFile.R);

            HashSet<String> tmpProbesList = new HashSet<String>(inProbes.readAsArrayList());

            TextFile e = new TextFile(inQTLs, TextFile.R);

            TextFile outWriter = new TextFile(out, TextFile.W);
            outWriter.writeln("Probe\tCombination\tP-value\tZscore\tFDR");

            String str = e.readLine();
            while((str = e.readLine())!=null){
                String[] parts = str.split(("\t"));

                if(tmpProbesList.contains(parts[4])){
                    outWriter.writeln(parts[4]+"\t"+parts[1]+"-"+parts[4]+"\t"+parts[0]+"\t"+parts[10]+"\t"+parts[(parts.length-1)]);
                }
            }

            outWriter.close();
        }
        
    }
}
