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

/**
 *
 * @author adriaan
 */
class CTSbinomialEntry {

    public CTSbinomialEntry(String asLocations,String phenoLocation, String outputLocation) throws FileNotFoundException, UnsupportedEncodingException, IOException, Exception {

        PrintWriter out_writer;
        out_writer = new PrintWriter(outputLocation, "UTF-8");

        // filenames_file, is a file containing per line a filename, which to load.
        String filenames_file = asLocations;

        /*
            Determine cell proportion per sample
         */
        
        ArrayList<String> phenoString = UtilityMethods.readFileIntoStringArrayList(phenoLocation);
        
        /*
            Right now just assuming this is a file that is 
            ordered in the same way as the asLocation file.
            With per line the cell proportion that we can determine.
            
            This is currently the users responsibility.
        */
        
        int i = 0;
        Double[] cellProp = new Double[phenoString.size()];
        for(String samplePheno : phenoString){
            cellProp[i] = Double.parseDouble(samplePheno);
            i++;
        }
        
        
        
        /*
            Open all the files we want to open, 
            read one line and combine the data 
            from the same SNP into the correct test.
        */
        
        ReadAsLinesIntoIndividualSNPdata asReader = new ReadAsLinesIntoIndividualSNPdata(asLocations);

        
        while (true) {
            
            //read some stuff from the files.
            ArrayList<IndividualSnpData> allSnpData;
            allSnpData = asReader.getIndividualsFromNextLine();
            
            if(allSnpData.isEmpty()) break;
            //Add the dispersion data assuming the same ordering
            //Which was checked previously.
            for(int j = 0; j < cellProp.length; j++){
                allSnpData.get(j).setCellTypeProp(cellProp[j]);
            }

            //do the binomial test:
            CTSbinomialTest results = new CTSbinomialTest(allSnpData);

            // Write the results to the out_file.
            if (results.isTestPerformed()) {
                out_writer.println(results.writeTestStatistics(true));
                GlobalVariables.numberOfTestPerformed++;
            }

        }

        out_writer.close();

        UtilityMethods.printFinalStats();
    }
    
}
