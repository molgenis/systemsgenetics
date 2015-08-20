/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;

/*
 *
 * @author Adriaan
 */

public class BinomialEntry {

    public static void BinomialEntry(String asLocations, String outputLocation) throws FileNotFoundException, UnsupportedEncodingException, IOException {

        PrintWriter out_writer;
        out_writer = new PrintWriter(outputLocation, "UTF-8");

        // filenames_file, is a file containing per line a filename, which to load.
        String filenames_file = asLocations;

        ReadAsLinesIntoIndividualSNPdata asReader = new ReadAsLinesIntoIndividualSNPdata(filenames_file);
        
        int i = 0;
        while(true){
            ArrayList<IndividualSnpData> allSnpData;
            //read all lines of the files.
            allSnpData = asReader.getIndividualsFromNextLine();
            
            if(allSnpData.isEmpty()) break; //when this occurs, all the files have no line
            
            
            BinomialTest results = new BinomialTest(allSnpData);
            
            // Write the results to the out_file.
            // And increase the number of tests.
            if (results.isTestPerformed()) {
                out_writer.println(results.writeTestStatistics(false));
                GlobalVariables.numberOfTestPerformed++;
            }
  
            
        }
        
        out_writer.close();

        UtilityMethods.printFinalTestStats();

    }




}
