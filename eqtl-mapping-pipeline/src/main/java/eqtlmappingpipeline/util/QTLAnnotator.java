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
        addAnnotationToQTLOutput(
                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Trans_Pc22c_CisMeQTLc_meQTLs\\RegressedOut_CisEffects\\eQTLsFDR0.05-ProbeLevel.txt",
                "D:\\UMCG\\GWAS_Catalog\\GWAS-Catalog-02092014.txt;D:\\UMCG\\Methylation_GPL13534\\annotationFile\\Illumina450K_MQtlMappingFile_Extensive.txt;D:\\UMCG\\Methylation_GPL13534\\annotationFile\\Illumina450K_MQtlMappingFile_Extensive.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\Optimal_PC_and_QTL_Corrected\\RP3_2MB_TSS_extendedCis_eQTMs\\eQTLSNPsFDR0.05-SNPLevel.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\Normal\\RP3_2Mb_extendedCis_eQTMs\\eQTLSNPsFDR0.05-SNPLevel.txt",
                "1;1;1;10;1;11;1;4;1;4", "snp;probe;probe;probe;probe", "D:\\UMCG\\ProbeMapping\\Info\\V70\\gencode.v15.annotation.gtf.gz",
                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Trans_Pc22c_CisMeQTLc_meQTLs\\RegressedOut_CisEffects\\test_eQTLsFDR0.05-ProbeLevel-ExtendedInfo.txt");

//        addAnnotationToQTLOutput(
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Cis_Pc22c_meQTLs\\Primary\\eQTLProbesFDR0.05-ProbeLevel_ldDrivenEffectsRemoved.txt",
//                "D:\\UMCG\\GWAS_Catalog\\GWAS-Catalog-02092014.txt;D:\\UMCG\\Methylation_GPL13534\\annotationFile\\Illumina450K_MQtlMappingFile_Extensive.txt;D:\\UMCG\\Methylation_GPL13534\\annotationFile\\Illumina450K_MQtlMappingFile_Extensive.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\Optimal_PC_and_QTL_Corrected\\RP3_2MB_TSS_extendedCis_eQTMs\\eQTLSNPsFDR0.05-SNPLevel.txt;D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\Normal\\RP3_2Mb_extendedCis_eQTMs\\eQTLSNPsFDR0.05-SNPLevel.txt", 
//                "1;1;1;10;1;11;1;4;1;4", "snp;probe;probe;probe;probe", null,
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Cis_Pc22c_meQTLs\\Primary\\eQTLProbesFDR0.05-ProbeLevel_ldDrivenEffectsRemoved-ExtendedInfo.txt");

//        addAnnotationToQTLOutput(
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\Optimal_PC_and_QTL_Corrected\\RP3_2MB_TSS_extendedCis_eQTMs\\eQTLSNPsFDR0.05-SNPLevel.txt",
//                "D:\\UMCG\\Methylation_GPL13534\\annotationFile\\Illumina450K_MQtlMappingFile_Extensive.txt", 
//                "1;10", "probe", null,
//                "D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\eQTMs\\Optimal_PC_and_QTL_Corrected\\RP3_2MB_TSS_extendedCis_eQTMs\\eQTLSNPsFDR0.05-SNPLevel.txt-ExtendedInfo.txt");
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
                extraHeaders += "\tAnnotionSource:" + masterId;
                try {
                    keyInt = Integer.parseInt(keyValueList[(masterId * 2)]);
                } catch (NumberFormatException ex) {
                    System.out.println("Error parsing key: " + keyValueList[masterId * 2] + " for annotationsource:" + annotationSource);
                    System.exit(-1);
                }

                int valueInt = 0;
                try {
                    valueInt = Integer.parseInt(keyValueList[((masterId * 2) + 1)]);
                } catch (NumberFormatException ex) {
                    System.out.println("Error parsing key: " + keyValueList[masterId * 2] + " for annotationsource:" + annotationSource);
                    System.exit(-1);
                }
                TextFile t = new TextFile(annotationSource, TextFile.R);
//                System.out.println(keyInt +"\t"+ valueInt);
                HashMap<String, String> annotationFromSource = (HashMap<String, String>) t.readAsHashMap(keyInt, valueInt);

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
//                        System.out.println("-");
                        if (!addTo) {
                            extraAnnotation.add("-");
                        } else {
                            extraAnnotation.set(id, extraAnnotation.get(id) + "\t-");
                        }

                    }
                    id++;

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
