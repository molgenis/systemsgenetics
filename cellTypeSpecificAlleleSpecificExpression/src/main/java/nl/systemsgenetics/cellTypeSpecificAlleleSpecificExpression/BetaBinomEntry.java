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
import org.apache.commons.io.FilenameUtils;
import org.jdom.IllegalDataException;

/**
 *
 * @author Adriaan van der Graaf 2015
 * 
 * This class will determine Maximum likelihood of Beta Binomial Over dispersion 
 * Then it will perform beta binomial likelihood ratio test per SNP.
 * 
 * Current uncertainty: 
 *                  I just take the likelihood of a certain a in null and compare 
 *                  this to freely estimated, don't know if this is correct.
 * 
 * Limitation to this is that it will read the as_files twice. 
 * Maybe there might be some better logic for this, but this would require quite 
 * some memory, so I'll read through them twice.
 * 
 */
public class BetaBinomEntry {
   
    
    public BetaBinomEntry(String asLocations, String outputLocation) throws FileNotFoundException, UnsupportedEncodingException, IOException, Exception{
        
        /*
         * PART 1: read all individual names
         */
        
        ArrayList<String> allFiles = UtilityMethods.readFileIntoStringArrayList(asLocations);
        
        
        /*
         * PART 2: Determine all Overdispersion parameters for all samples.
         * TODO Implement some logic if there is already overdispersion data available.
         *      Perhaps implement this inside the BetaBinomOverdispInSample side.
         *      Not necassary now, but later, yes.
         */
        
        ArrayList<BetaBinomOverdispInSample>  dispersionParameters = new ArrayList<BetaBinomOverdispInSample>();
        
        for(String sampleName : allFiles){
           
            dispersionParameters.add(new BetaBinomOverdispInSample(sampleName));
            /*
             * TODO: determine some way to save the overdispersion sample.
             *       But we could also do it everytime you do this, perhaps based 
             *       on the relative time compared to the full calculation.
             *       
                    EDIT: seems to be not too long (less than a minute) per individual, 
             *       so I'm fine with keeping it this way.
             */
        }
        double[] dispersionVals = new double[dispersionParameters.size()];  
        int i = 0;     
        
        
        // use this to save the dispersionvalues
        
        String dispersionOutput = FilenameUtils.getPath(outputLocation) + 
                                  FilenameUtils.getBaseName(outputLocation) +
                                  "_dispersionFile.txt";
        
        PrintWriter writer = new PrintWriter(dispersionOutput, "UTF-8");       
        
        //header
        writer.write("Filename\tdispersion");
        
        
        
        for(BetaBinomOverdispInSample sampleDispersion : dispersionParameters){
        
            dispersionVals[i] = sampleDispersion.getOverdispersion()[0];
            
            //do a check to make sure ordering is correct.
            if(!(sampleDispersion.getSampleName().equals(allFiles.get(i)))){
                System.out.println(sampleDispersion.getSampleName());
                throw new IllegalDataException("ERROR! ordering is not correct filenames for overdispersion");
            }
            
            writer.printf("%s\t%.6f\n", sampleDispersion.getSampleName(), sampleDispersion.getOverdispersion()[0] );
            
            i++;
        }
        
        writer.close();
        
        
        

                
                
        /*
         * Part 3: Read the AS_files per line, and do the Binomial Test per SNP.
         *         
         */
        
        
        PrintWriter out_writer;
        out_writer = new PrintWriter(outputLocation, "UTF-8");
        
        //open all the files we want to open.

        ReadAsLinesIntoIndividualSNPdata asReader = new ReadAsLinesIntoIndividualSNPdata(asLocations);
        
        while (true) {
            //read some stuff from the files.
            ArrayList<IndividualSnpData> allSnpData;
            allSnpData = asReader.getIndividualsFromNextLine();
            
            if(allSnpData.isEmpty()) break;
            //Add the dispersion data assuming the same ordering
            //Which was checked previously.
            for(int j = 0; j < dispersionVals.length; j++){
                allSnpData.get(j).setDispersion(dispersionVals[j]);
            }

            
            //do the beta binomial test:
            BetaBinomialTest results = new BetaBinomialTest(allSnpData);

            // Write the results to the out_file.
            if (results.isTestPerformed()) {
                out_writer.println(results.writeTestStatistics(true));
                GlobalVariables.numberOfTestPerformed++;
            }

        }
        
        out_writer.close();
        UtilityMethods.printFinalTestStats();

        
    }

}
