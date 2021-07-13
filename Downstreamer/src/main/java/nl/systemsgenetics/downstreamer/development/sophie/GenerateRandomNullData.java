/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.development.sophie;

/**
 *
 * @author Sophie Mulc
 *
 */
import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import com.opencsv.CSVWriter;
import java.io.BufferedReader;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

@Deprecated
public class GenerateRandomNullData {

    /**
     * @param args the command line arguments
     * @throws java.io.FileNotFoundException
     */
    public static void main(String[] args) throws FileNotFoundException, IOException {
        // TODO code application logic here

        File phase3File = new File("C:\\Users\\Sophie Mulc\\Documents\\DEPICT2\\phase3_corrected.psam");
        File phase3RandomFile = new File("C:\\Users\\Sophie Mulc\\Documents\\DEPICT2\\phase3_corrected_randomPhenosNew.psam");

        //FileReader(String phase3_corrected)
        final CSVParser gmtParser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
        final CSVReader gmtReader = new CSVReaderBuilder(new BufferedReader(new FileReader(phase3File))).withSkipLines(1).withCSVParser(gmtParser).build();

        CSVWriter writer = new CSVWriter(new FileWriter(phase3RandomFile), '\t', '\0', '\0', "\n");

     
        
        String[] outputLine = new String[1002];
        int c = 0;
        outputLine[c++] = "#FID";
        outputLine[c++] = "#IID";
        for (int i = 2; i < 1002; ++i) {
        outputLine[c++] = "RA" + i;
        }
        writer.writeNext(outputLine);

        List<String> iids = new ArrayList<>();

        String[] nextLine;
        while ((nextLine = gmtReader.readNext()) != null) {

            String iid = nextLine[0];

            iids.add(iid);

        }

        System.out.println(iids);

        System.out.println("IIDS");
        
          Random randomno1 = new Random();
          
        for (String iid : iids) {
            System.out.println(iid);

            c = 0;
            outputLine[c++] = iid;
            outputLine[c++] = iid;
            
            for (int i = 2; i < 1002; ++i) {
            outputLine[c++] = String.valueOf(randomno1.nextGaussian());
            }
            writer.writeNext(outputLine);
            

        }
    
      writer.close();
//   
    }

}
