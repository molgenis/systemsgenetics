/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.util;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author MarcJan
 */
public class SummarizeQtlSummaries {
    public static void main(String[] args) {
        File basePath = new File("D:\\OnlineFolders\\AeroFS\\SharedQTLFolder\\QTL_Output\\For_paper_20160208\\Replication");
        String filenameContains= "-Summary.txt";
        
        ArrayList<File> fileList = new ArrayList<>();
        listFilesForFolder(basePath, fileList);
        
        try {
            processData(fileList, filenameContains);
        } catch (IOException ex) {
            Logger.getLogger(SummarizeQtlSummaries.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    }
    
    public static void listFilesForFolder(File folder, ArrayList<File> fileList) {
        for (final File fileEntry : folder.listFiles()) {
            if (fileEntry.isDirectory()) {
                listFilesForFolder(fileEntry, fileList);
            } else {
                fileList.add(fileEntry);
            }
        }
    }

    private static void processData(ArrayList<File> fileList, String filenameContains) throws IOException {
        for(File f : fileList){
            if(f.getAbsolutePath().endsWith(filenameContains)){
                String fileName = f.getAbsolutePath();
                TextFile t = new TextFile(fileName, TextFile.R);
//                System.out.println(fileName);
                String str;
                
                int counter = 0;
                int eQTLsInDisc = 0;
                String file1 = "";
                int eQTLsInRep = 0;
                String file2 = "";
                int qtlsShared = 0; 
                int nrIdentical = 0;
                int nrOpposite = 0;
                while((str=t.readLine())!=null){
                    String[] parts = str.split("\t");

                    if(counter ==0){
                        eQTLsInDisc = Integer.parseInt(parts[1]);
                        file1 = parts[3];
                        file1 = file1.replace("D:\\OnlineFolders\\AeroFS\\SharedQTLFolder\\QTL_Output\\For_paper_20160208\\LifeLines_Discovery\\","");
                    } else if(counter == 1){
                        eQTLsInRep = Integer.parseInt(parts[1]);
                        file2 = parts[3];
                        file2 = file2.replace("D:\\OnlineFolders\\AeroFS\\SharedQTLFolder\\QTL_Output\\For_paper_20160208\\Replication\\","");
                    } else if(counter == 2){
                        qtlsShared = Integer.parseInt(parts[1]);
                    } else if(counter == 4){
                        nrIdentical = Integer.parseInt(parts[1]);
                    }  else if(counter == 5){
                        nrOpposite = Integer.parseInt(parts[1]);
                    }  

//                    System.out.println("\t"+str);
                    counter++;
                }
                System.out.println(f.getName().replace("-Summary.txt", "")+"\t"+file1+"\t"+eQTLsInDisc+"\t"+file2+"\t"+eQTLsInRep+"\t"+qtlsShared+"\t"+nrIdentical+"\t"+nrOpposite);
            }
        }
    }
}
