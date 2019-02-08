/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2.development;

import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import java.io.BufferedReader;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;
import static org.apache.hadoop.hive.common.StatsSetupConst.StatDB.counter;

/**
 *
 * @author Sophie Mulc
 */
public class BinaryCompareFormat {

    /**
     * @param args the command line arguments
     * @throws java.io.IOException
     */
    public static void main(String[] args) throws IOException {
        // TODO code application logic here

        //To compare the format of the files we created in the eQTL pipeline for running GWAS to the format needed to run in Lude's code.
        final File predictionMatrixFile = new File("/groups/umcg-wijmenga/scr02/umcg-smulcahy/eQTLResults_texttrue2/eQTL.binary");
        System.out.println(predictionMatrixFile);
        DoubleMatrixDataset<String, String> predictionMatrixFull = DoubleMatrixDataset.loadTransEqtlExpressionMatrix(predictionMatrixFile.getAbsolutePath());

        System.out.println(predictionMatrixFull.getElement("rs351365", "PH443"));
        File eQTLfile = new File("/groups/umcg-wijmenga/scr02/umcg-smulcahy/eQTLResults_texttrue2/eQTLs.txt.gz");
        final CSVParser eQTLparser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
        final CSVReader eQTLreader = new CSVReaderBuilder((new InputStreamReader(new GZIPInputStream(new FileInputStream(eQTLfile))))).withSkipLines(1).withCSVParser(eQTLparser).build();

        int c = 0;
        String[] inputLine;
        while ((inputLine = eQTLreader.readNext()) != null) {

            String snp = inputLine[1];
            String pheno = inputLine[4];
            String zscore = inputLine[10];

            // Test with one site only need this line and to read in file: predictionMatrixFull.getElement("rs351365", "PH443"); i.e (rowName(snp), columnName(pheno))
        

            double zscore2 = (predictionMatrixFull.getElement(snp, pheno));

            //convert string to double, then look up how to compare two doubles - this is to compare the zscores
            double d_zscore = Double.parseDouble(zscore);

            double compare_zscores = d_zscore - zscore2;

            //count occurances of above 0 comparisons
            if (Math.abs(compare_zscores) > 0.00001) {
                c++;
            }
        }
        System.out.println("Number of occurrances where z-scores differ: " + c);
        //save new file 
        predictionMatrixFull.saveBinary("/groups/umcg-wijmenga/scr02/umcg-smulcahy/eQTLResults_texttrue2/eQTL2");
        
        System.out.println("Done saving");

    }

}
