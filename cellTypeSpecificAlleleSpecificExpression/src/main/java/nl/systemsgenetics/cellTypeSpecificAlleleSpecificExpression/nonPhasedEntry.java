/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import org.apache.commons.io.FilenameUtils;
import org.jdom.IllegalDataException;

/**
 *
 * @author adriaan
 */
public class nonPhasedEntry {
    
    
    public nonPhasedEntry(String asLocations, String phenoTypeLocation, String outputLocation) throws IOException, Exception {
    
        //This is the entry used for non-phased tests.
        
        //PART 1: read all individuals names from the files.

        ArrayList<String> allFiles = UtilityMethods.readFileIntoStringArrayList(asLocations);
        
        //PART 2: determine the per sample overdispersion in the file.
        
        ArrayList<BetaBinomOverdispInSample>  dispersionParameters = new ArrayList<BetaBinomOverdispInSample>();

        
        
        
        for(String sampleName : allFiles){
            dispersionParameters.add(new BetaBinomOverdispInSample(sampleName));
        }
        

        // use this to save the dispersionvalues
        
        String dispersionOutput = FilenameUtils.getPath(outputLocation) + 
                                  FilenameUtils.getBaseName(outputLocation) +
                                  "_dispersionFile.txt";
        
        String binomialOutput  = FilenameUtils.getPath(outputLocation) + 
                                  FilenameUtils.getBaseName(outputLocation) +
                                  "_BinomialResults.txt";
        
        String betaBinomialOutput = FilenameUtils.getPath(outputLocation) + 
                                    FilenameUtils.getBaseName(outputLocation) +
                                    "_BetaBinomialResults.txt";
        
        //for ease of use initializing here.
        String CTSbinomialOutput  = FilenameUtils.getPath(outputLocation) + 
                                  FilenameUtils.getBaseName(outputLocation) +
                                  "_CTSBinomialResults.txt";
        
        String CTSbetaBinomialOutput = FilenameUtils.getPath(outputLocation) + 
                                    FilenameUtils.getBaseName(outputLocation) +
                                    "_CTSBetaBinomialResults.txt";
        
        
        
        
        
        PrintWriter writer = new PrintWriter(dispersionOutput, "UTF-8");       
        
        //header for dispersion
        writer.write("Filename\tdispersion");        
        

        double[] dispersionVals = new double[dispersionParameters.size()];  
        int i = 0;     
        
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
        
        
        //This is  only done when there is a correct location of the pheno file
        double[] cellProp = new double[] {-1.0};
        
        if(phenoTypeLocation != null){
            ArrayList<String> phenoString =UtilityMethods.readFileIntoStringArrayList(phenoTypeLocation);

            /*
                Right now just assuming this is a file that is 
                ordered in the same way as the asLocation file.
                With per line the cell proportion that we can determine.
                This is a requirement of the input!
            */
            i = 0;
            cellProp = new double[phenoString.size()];
            for(String samplePheno : phenoString){
                cellProp[i] = Double.parseDouble(samplePheno);
                i++;
            }
        }
        
        //PART 4. Read the as files one line at a time.
        
        //Will create three types of 
        
        PrintWriter binomWriter = new PrintWriter(binomialOutput, "UTF8");
        PrintWriter betaBinomWriter = new PrintWriter(betaBinomialOutput, "UTF8");

        
        
        //CTS stuff.
        PrintWriter CTSBinomWriter = null; 
        PrintWriter CTSBetaBinomWriter = null;
        
        if(phenoTypeLocation != null){
            CTSBinomWriter = new PrintWriter(CTSbinomialOutput, "UTF8");
            CTSBetaBinomWriter = new PrintWriter(CTSbetaBinomialOutput, "UTF-8");
        }
        
        
        //open all the files we want to open.
        ReadAsLinesIntoIndividualSNPdata asReader = new ReadAsLinesIntoIndividualSNPdata(asLocations);

        while (true) {

            //read some the next line from the files.
            ArrayList<IndividualSnpData> allSnpData;
            allSnpData = asReader.getIndividualsFromNextLine();
            
            
            //BREAKPOINT OF THE LOOP!!!
            if(allSnpData.isEmpty()) break;
            
            
            
            // Add the dispersion data assuming the same ordering
            // Which was checked previously.
            for(int j = 0; j < dispersionVals.length; j++){
                allSnpData.get(j).setDispersion(dispersionVals[j]);
            }
            
            
            BinomialTest binomialResults = new BinomialTest(allSnpData);
            BetaBinomialTest betaBinomialResults = new BetaBinomialTest(allSnpData);
            
            if(binomialResults.isTestPerformed()){
            
                binomWriter.println(binomialResults.writeTestStatistics(false));
                betaBinomWriter.println(betaBinomialResults.writeTestStatistics(false));
                
                GlobalVariables.numberOfTestPerformed++;
                
            }
            
            

            // add the cellProp to the snp and do the CTS tests.
            if(phenoTypeLocation != null){
                for(int j = 0; j < cellProp.length; j++){
                    allSnpData.get(j).setCellTypeProp(cellProp[j]);
                }
                //do the CTS beta binomial test:

                CTSBetaBinomialTest CTSbetaBinomResults = new CTSBetaBinomialTest(allSnpData);
                CTSbinomialTest CTSbinomialResults = new CTSbinomialTest(allSnpData);
                
                // Write the results to the out_file, assuming both of them were done.
                if (CTSbetaBinomResults.isTestPerformed()) {
                    
                    CTSBinomWriter.println(CTSbinomialResults.writeTestStatistics(true));
                    CTSBetaBinomWriter.println(CTSbetaBinomResults.writeTestStatistics(true));
                }
            }
        
        }
     
        binomWriter.close();
        betaBinomWriter.close();
        
        
        if(phenoTypeLocation != null){    
            CTSBinomWriter.close();
            CTSBetaBinomWriter.close();
        }
        
        UtilityMethods.printFinalTestStats();
        
    }


}
