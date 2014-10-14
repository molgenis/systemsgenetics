/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.util;

import java.io.IOException;
import java.util.ArrayList;
import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;

/**
 *
 * @author MarcJan
 */
class QTLAnnotator {
    
    public static void main(String[] args) throws IOException {
        addAnnotationToQTLOutput("D:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Trans_Pc22c_CisMeQTLc_meQTLs\\RegressedOut_CisEffects\\eQTLsFDR0.05-ProbeLevel.txt", 
                "D:\\UMCG\\GWAS_Catalog\\GWAS-Catalog-02092014.txt", "1;1");
    }

    static void addAnnotationToQTLOutput(String in, String sources, String keyValuePairs) throws IOException {
        QTLTextFile e = new QTLTextFile(in, QTLTextFile.R);
        
        ArrayList<EQTL> qtls = e.readList();
        String extraHeaders = "";
        ArrayList<String> extraAnnotation = new ArrayList<String>(qtls.size());
        
        String[] sourceList = sources.split(";");
        String[] keyValueList = keyValuePairs.split(";");
        
        for(String annotationSource : sourceList){
            if(annotationSource.contains("GWAS") &&annotationSource.contains("Catalog")){
                GWASCatalog cat = new GWASCatalog(annotationSource);
                
                for(EQTL eQtl : qtls){
                    System.out.println(cat.getFullTraitsForCertainSnps(eQtl.getRsName()));
                }
                
            } else {
                
            }
            
        }
        
    }
    
}
