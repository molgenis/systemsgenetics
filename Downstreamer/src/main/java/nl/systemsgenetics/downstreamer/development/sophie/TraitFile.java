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
public class TraitFile {

    /**
     * @param args the command line arguments
     * @throws java.io.FileNotFoundException
     */
    public static void main(String[] args) throws FileNotFoundException, IOException {
        // TODO code application logic here

//        MAKE TRAIT FILE
        File phase3File = new File("C:\\Users\\Sophie Mulc\\Documents\\DEPICT2\\phase3_corrected.psam");
        File traitFile = new File("C:\\Users\\Sophie Mulc\\Documents\\DEPICT2\\TraitFile.txt");

        //FileReader(String phase3_corrected)
        final CSVParser gmtParser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
        final CSVReader gmtReader = new CSVReaderBuilder(new BufferedReader(new FileReader(phase3File))).withSkipLines(1).withCSVParser(gmtParser).build();

        List<String> iids = new ArrayList<>();

        String[] inputLine;
        while ((inputLine = gmtReader.readNext()) != null) {

            String iid = inputLine[0];

            iids.add(iid);

        }

        CSVWriter writer = new CSVWriter(new FileWriter(traitFile), '\t', '\0', '\0', "\n");

        String[] outLine = new String[iids.size() + 1];
        int c = 0;
        outLine[c++] = "PHENO";
        for (String iid : iids) {
            outLine[c++] = iid;
        }
        writer.writeNext(outLine);

        System.out.println(iids);

        Random randomno1 = new Random();

        for (int j = 1; j < 1001; ++j) {
            c = 0;

            outLine[c++] = "PH" + j;
            for (int i = 0; i < iids.size(); ++i) {
                outLine[c++] = String.valueOf(randomno1.nextGaussian());
            }
            writer.writeNext(outLine);

        }
        writer.close();

//  MAKE PROBE ANNOTATION FILE      
        File probeAnnotationFile = new File("C:\\Users\\Sophie Mulc\\Documents\\DEPICT2\\ProbeAnnotationFile.txt");

        CSVWriter writer1 = new CSVWriter(new FileWriter(probeAnnotationFile), '\t', '\0', '\0', "\n");

        String[] output1Line = new String[6];
        c = 0;
        output1Line[c++] = "Platform";
        output1Line[c++] = "HT12v4-ArrayAddress";
        output1Line[c++] = "Symbol";
        output1Line[c++] = "Chr";
        output1Line[c++] = "ChrStart";
        output1Line[c++] = "ChrEnd";
        writer1.writeNext(output1Line);

        for (int j = 1; j < 1001; ++j) {
            c = 0;

            output1Line[c++] = "Plat";
            output1Line[c++] = "PH" + j;
            output1Line[c++] = "PH" + j;
            output1Line[c++] = "25";
            output1Line[c++] = "100";
            output1Line[c++] = "154";
            writer1.writeNext(output1Line);
        }
        writer1.close();

        //  MAKE COUPLING FILE      
        File couplingFile = new File("C:\\Users\\Sophie Mulc\\Documents\\DEPICT2\\CouplingFile.txt");

        CSVWriter writer2 = new CSVWriter(new FileWriter(couplingFile), '\t', '\0', '\0', "\n");

        String[] output2Line = new String[2];
        for (String iid : iids) {
            c = 0;
            output2Line[c++] = iid;
            output2Line[c++] = iid;
            writer2.writeNext(output2Line);
        }

        writer2.close();
    }

}
