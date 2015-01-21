/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.util;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.io.gwascatalog.GWASSNP;
import umcg.genetica.io.gwascatalog.GWASTrait;
import umcg.genetica.io.probemapping.reading;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;

/**
 *
 * @author MarcJan
 */
class QTLAnnotator {

    public static void main(String[] args) throws IOException {
//        addAnnotationToQTLOutput(
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Trans_Pc22c_CisMeQTLc_meQTLs\\RegressedOut_CisEffects\\eQTLsFDR0.05-ProbeLevel.txt",
//                "D:\\UMCG\\GWAS_Catalog\\GWAS-Catalog-02092014.txt;D:\\UMCG\\Methylation_GPL13534\\annotationFile\\Illumina450K_MQtlMappingFile_Extensive.txt;D:\\UMCG\\Methylation_GPL13534\\annotationFile\\Illumina450K_MQtlMappingFile_Extensive.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\QTLCorrected\\RP3_1MB_TSS_extendedCis_eQTMs\\eQTLSNPsFDR0.05-Closest-SNPLevel.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\QTLCorrected\\RP3_2MB_TSS_extendedCis_eQTMs\\eQTLSNPsFDR0.05-Closest-SNPLevel.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\Great_GeneAssociations_20141127-public-2.0.2-aaQ50i-hg19-all-region_Ensembl.txt",
//                "1;1;1;10;1;11;1;4;1;4;0;1", "snp;probe;probe;probe;probe;probe", "D:\\UMCG\\ProbeMapping\\Info\\V70\\gencode.v15.annotation.gtf.gz",
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Trans_Pc22c_CisMeQTLc_meQTLs\\RegressedOut_CisEffects\\eQTLsFDR0.05-ProbeLevel-ExtendedInfoTesting.txt");
        
//        addAnnotationToQTLOutput(
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Trans_Pc22c_CisMeQTLc_meQTLs\\RegressedOut_CisEffects_New\\eQTLs_head.txt",
//                "D:\\UMCG\\GWAS_Catalog\\GWAS-Catalog-02092014.txt;D:\\UMCG\\Methylation_GPL13534\\annotationFile\\Illumina450K_MQtlMappingFile_Extensive.txt;D:\\UMCG\\Methylation_GPL13534\\annotationFile\\Illumina450K_MQtlMappingFile_Extensive.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\QTLCorrected\\RP3_1MB_TSS_extendedCis_eQTMs\\eQTLSNPsFDR0.05-Closest-SNPLevel.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\QTLCorrected\\RP3_2MB_TSS_extendedCis_eQTMs\\eQTLSNPsFDR0.05-ClosestStrongest-SNPLevel.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\Great_GeneAssociations_20141127-public-2.0.2-aaQ50i-hg19-all-region_Ensembl.txt",
//                "1;1;1;10;1;11;1;4;1;4;0;1", "snp;probe;probe;probe;probe;probe", "D:\\UMCG\\ProbeMapping\\Info\\V70\\gencode.v15.annotation.gtf.gz",
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Trans_Pc22c_CisMeQTLc_meQTLs\\RegressedOut_CisEffects_New\\eQTLs_head-ExtendedInfoTesting.txt");

//        addAnnotationToQTLOutput(
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Cis_Pc22c_meQTLs\\Primary\\eQTLProbesFDR0.05-ProbeLevel_ldDrivenEffectsRemoved.txt",
//                "D:\\UMCG\\GWAS_Catalog\\GWAS-Catalog-02092014.txt;D:\\UMCG\\Methylation_GPL13534\\annotationFile\\Illumina450K_MQtlMappingFile_Extensive.txt;D:\\UMCG\\Methylation_GPL13534\\annotationFile\\Illumina450K_MQtlMappingFile_Extensive.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\Optimal_PC_and_QTL_Corrected\\RP3_2MB_TSS_extendedCis_eQTMs\\eQTLSNPsFDR0.05-SNPLevel.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\Normal\\RP3_2Mb_TSSextendedCis_eQTMs\\eQTLSNPsFDR0.05-SNPLevel.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\QTLCorrected\\RP3_2MB_TSS_extendedCis_eQTMs\\eQTLSNPsFDR0.05-SNPLevel.txt",
//                "1;1;1;10;1;11;1;4;1;4;1;4", "snp;probe;probe;probe;probe;probe", "D:\\UMCG\\ProbeMapping\\Info\\V70\\gencode.v15.annotation.gtf.gz",
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Cis_Pc22c_meQTLs\\Primary\\eQTLProbesFDR0.05-ProbeLevel_ldDrivenEffectsRemoved-ExtendedInfo.txt");

//        addAnnotationToQTLOutput(
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\Optimal_PC_and_QTL_Corrected\\RP3_0.25MB_TSS_extendedCis_eQTMs_2015\\eQTLSNPsFDR0.05-SNPLevel.txt",
//                "D:\\UMCG\\Methylation_GPL13534\\annotationFile\\Illumina450K_MQtlMappingFile_Extensive.txt;D:\\UMCG\\Data\\RP3_RNA_Seq\\annotation_geneIds+overlapping_v71_cut.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\Annotation450k.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\CODAM_NTR_LLS_LLD_RS_BBMRI_450K_var_mean_median.txt", 
//                "1;8-9-10-11-12-13-14;1;4-5;0;17-18-20-49-50-51-52-53-54-55-56-57-58-59-60-61-62-63-64-65-66-67-68-70-71-81-133-134-141-142;0;1-2-3-4", "snp;probe;snp;snp", null,
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\Optimal_PC_and_QTL_Corrected\\RP3_0.25MB_TSS_extendedCis_eQTMs_2015\\eQTLSNPsFDR0.05-SNPLevel.txt-ExtendedInfo2.txt");
        
        addAnnotationToQTLOutput(
                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\QTLCorrected\\RP3_0.25MB_TSS_extendedCis_eQTMs_2015\\eQTLSNPsFDR0.05-SNPLevel.txt",
                "D:\\UMCG\\Methylation_GPL13534\\annotationFile\\Illumina450K_MQtlMappingFile_Extensive.txt;D:\\UMCG\\Data\\RP3_RNA_Seq\\annotation_geneIds+overlapping_v71_cut.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\Annotation450k_AdditionMJ.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\CODAM_NTR_LLS_LLD_RS_BBMRI_450K_var_mean_median.txt", 
                "1;8-9-10-11-12-13-14;1;4-5;0;17-18-20-49-50-51-52-53-54-55-56-57-58-59-60-61-62-63-64-65-66-67-68-70-71-81-133-134-141-142-148;0;1-2-3-4", "snp;probe;snp;snp", null,
                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\QTLCorrected\\RP3_0.25MB_TSS_extendedCis_eQTMs_2015\\eQTLSNPsFDR0.05-SNPLevel_ExtendedInfo3.txt");
    }

    static void addAnnotationToQTLOutput(String in, String sources, String keyValuePairs, String idsToAnnotate, String reannotateGene, String out) throws IOException {
        QTLTextFile e = new QTLTextFile(in, QTLTextFile.R);

        ArrayList<EQTL> qtls = e.readList();
        String extraHeaders = "";
        ArrayList<String> extraAnnotation = new ArrayList<String>(qtls.size());

        String[] sourceList = sources.split(";");
        String[] keyValueList = keyValuePairs.split(";");
        String[] annotationIdList = idsToAnnotate.split(";");

        if (sourceList.length * 2 != keyValueList.length) {
            System.out.println("Faulty key value match input.");
            System.exit(0);
        }
        
        
        boolean addTo = false;
        
        if (!(reannotateGene == null || reannotateGene.equals(""))) {
            extraHeaders += "\tEnsemble Gene ids (;)";
            HashMap<String, String> gencodeGeneMapping = reading.readGTFAnnotationFileHash(reannotateGene, 50000);
            HashMap<String, HashSet<String>> genCodeMappingGenes = new HashMap<String, HashSet<String>>();
            for (Map.Entry<String, String> t : gencodeGeneMapping.entrySet()) {
                if (t.getKey().contains("ENSG")) {
                    String[] newKeys;
                    if (t.getValue().contains(";")) {
                        newKeys = t.getValue().split(";");
                    } else {
                        newKeys = new String[1];
                        newKeys[0] = t.getValue();
                    }

                    for (String nKey : newKeys) {
                        if (genCodeMappingGenes.containsKey(nKey)) {
                            genCodeMappingGenes.get(nKey).add(t.getKey());
                        } else {
                            HashSet<String> ensembleIds = new HashSet<String>();
                            ensembleIds.add(t.getKey());
                            genCodeMappingGenes.put(nKey, ensembleIds);
                        }
                    }
                }
            }
            int id = 0;
            for (EQTL qtl : qtls) {
                qtl.setProbeHUGO(qtl.getProbeHUGO().replace("\"", ""));
                
                if (qtl.getProbeHUGO().contains(";")) {
                    String[] ts = qtl.getProbeHUGO().split(";");
                    HashSet<String> setje = new HashSet<String>();
                    for (String x : ts) {
                        if (genCodeMappingGenes.containsKey(x)) {
                            setje.addAll(genCodeMappingGenes.get(x));
                        }
                    }
                    if(setje.isEmpty()){
                        if (!addTo) {
                            extraAnnotation.add("-");
                        } else {
                            extraAnnotation.set(id, extraAnnotation.get(id) + "\t-" );
                        }
                    } else {
                        if (!addTo) {
                            extraAnnotation.add(setje.toString().replaceAll("\\[", "").replaceAll("\\]", "").replaceAll(", ", ";"));
                        } else {
                            extraAnnotation.set(id, extraAnnotation.get(id) + "\t" + setje.toString().replaceAll("\\[", "").replaceAll("\\]", "").replaceAll(", ", ";"));
                        }
                    }
                    
                } else {
                    if (genCodeMappingGenes.containsKey(qtl.getProbeHUGO())) {
                        if (!addTo) {
                            extraAnnotation.add(genCodeMappingGenes.get(qtl.getProbeHUGO()).toString().replaceAll("\\[", "").replaceAll("\\]", "").replaceAll(", ", ";"));
                        } else {
                            extraAnnotation.set(id, extraAnnotation.get(id) + "\t" + genCodeMappingGenes.get(qtl.getProbeHUGO()).toString().replaceAll("\\[", "").replaceAll("\\]", "").replaceAll(", ", ";"));
                        }
                    } else {
                        if (!addTo) {
                            extraAnnotation.add("-");
                        } else {
                            extraAnnotation.set(id, extraAnnotation.get(id) + "\t-" );
                        }
                    }
                }
                id++;
            }
        }

        int masterId = 0;
        for (String annotationSource : sourceList) {
            if (!extraAnnotation.isEmpty()) {
                addTo = true;
            }
            if (annotationSource.contains("GWAS") && annotationSource.contains("Catalog")) {

                GWASCatalog cat = new GWASCatalog(annotationSource);
                HashMap<String, GWASSNP> catSnps = cat.getSnpToObj();
                extraHeaders += "\tGWASTraits (;)\tRisk allels (;)\t P-value GWAS (;)";

                int id = 0;
                for (EQTL eQtl : qtls) {
                    String gwasTrait = "";
                    String gwasAllelle = "";
                    String gwasPval = "";
                    GWASSNP s = catSnps.get(eQtl.getRsName());
                    if (s != null) {
                        HashSet<GWASTrait> traits = s.getAssociatedTraits();

                        for (GWASTrait t : traits) {
                            if (!gwasTrait.equals("")) {
                                gwasTrait += ";";
                                gwasAllelle += ";";
                                gwasPval += ";";
                            }
                            gwasTrait += t.getCleanName();
                            gwasAllelle += s.getRiskAllele(t);
                            gwasPval += s.getPValueAssociatedWithTrait(t);

                        }
//                        System.out.println(gwasTrait + "\t" + gwasAllelle+ "\t" + gwasPval);
                        if (!addTo) {
                            extraAnnotation.add(gwasTrait + "\t" + gwasAllelle + "\t" + gwasPval);
                        } else {
                            extraAnnotation.set(id, extraAnnotation.get(id) + "\t" + gwasTrait + "\t" + gwasAllelle + "\t" + gwasPval);
                        }


                    } else {
//                        System.out.println("-\t-\t-");
                        if (!addTo) {
                            extraAnnotation.add("-\t-\t-");
                        } else {
                            extraAnnotation.set(id, extraAnnotation.get(id) + "\t-\t-\t-");
                        }

                    }
                    id++;
                }

            } else {
//                System.out.println(annotationSource);
//                System.out.println(masterId);
                int keyInt = 0;
                
                try {
                    keyInt = Integer.parseInt(keyValueList[(masterId * 2)]);
                } catch (NumberFormatException ex) {
                    System.out.println("Error parsing key: " + keyValueList[masterId * 2] + " for annotationsource:" + annotationSource);
                    System.exit(-1);
                }

                int[] valueInt = {0};
                try {
                    if(keyValueList[((masterId * 2) + 1)].contains("-")){
                        String[] tmpValues = keyValueList[((masterId * 2) + 1)].split("-");
                        valueInt = new int[tmpValues.length];
                        for(int i=0; i<valueInt.length; ++i){
                            valueInt[i] = Integer.parseInt(tmpValues[i]);
                        }
                    } 
                    else {
                        valueInt[0] = Integer.parseInt(keyValueList[((masterId * 2) + 1)]);
                    }
                    
                } catch (NumberFormatException ex) {
                    System.out.println("Error parsing key: " + keyValueList[masterId * 2] + " for annotationsource:" + annotationSource);
                    System.exit(-1);
                }
                
//                System.out.println(keyInt +"\t"+ valueInt);
                for(int valInt : valueInt){
                    extraHeaders += "\tAnnotionSource:" + masterId+"_Col:"+valInt;
                    TextFile t = new TextFile(annotationSource, TextFile.R);
                    HashMap<String, String> annotationFromSource = (HashMap<String, String>) t.readAsHashMap(keyInt, valInt);

                    int id = 0;
                    for (EQTL eQtl : qtls) {

                        String newInfo;
                        if (annotationIdList[masterId].equalsIgnoreCase("probe")) {
                            newInfo = annotationFromSource.get(eQtl.getProbe());
                        } else {
                            newInfo = annotationFromSource.get(eQtl.getRsName());
                        }

                        if (newInfo != null) {
    //                        System.out.println(newInfo);
                            if (!addTo) {
                                extraAnnotation.add(newInfo);
                            } else {
                                extraAnnotation.set(id, extraAnnotation.get(id) + "\t" + newInfo);
                            }

                        } else {
                            System.out.println("here");
                            System.exit(0);
    //                        System.out.println("-");
                            if (!addTo) {
                                extraAnnotation.add("-");
                            } else {
                                
                                extraAnnotation.set(id, extraAnnotation.get(id) + "\t-");
                            }

                        }
                        id++;
                    }
                    if (!extraAnnotation.isEmpty()) {
                        addTo = true;
                    }
                    t.close();
                }
            }
            masterId++;
        }

        TextFile outWriter = new TextFile(out, TextFile.W);

        outWriter.writeln(QTLTextFile.header + extraHeaders);
        int id = 0;
        for (EQTL qtl : qtls) {
            outWriter.writeln(qtl.toString() + "\t" + extraAnnotation.get(id));
            id++;
        }
        outWriter.close();
    }
}
