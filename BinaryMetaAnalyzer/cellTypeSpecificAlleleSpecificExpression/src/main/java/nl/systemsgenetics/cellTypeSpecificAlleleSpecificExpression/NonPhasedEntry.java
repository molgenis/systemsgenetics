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
public class NonPhasedEntry {
    
    
    public NonPhasedEntry(String asLocations, String phenoTypeLocation, String dispersionLocation, String outputLocation) throws IOException, Exception {
    
        //This is the entry used for non-phased tests.
        
        //PART 1: read all individuals names from the files.

        ArrayList<String> allFiles = UtilityMethods.readFileIntoStringArrayList(asLocations);
        
        //PART 2: determine the per sample overdispersion in the file.
        
        ArrayList<BetaBinomOverdispInSample>  dispersionParameters = new ArrayList<BetaBinomOverdispInSample>();

        
        
        if(dispersionLocation == null){
            for(String sampleName : allFiles){
                dispersionParameters.add(new BetaBinomOverdispInSample(sampleName));
            }
        }else{
            ArrayList<String> DispersionLines = UtilityMethods.readFileIntoStringArrayList(dispersionLocation);
            //remove first Dispersion
            DispersionLines.remove(0);
            for(String Line : DispersionLines){
                String sampleName = Line.split("\t")[0];
                double[] dispersion;
                dispersion = new double[] { Double.parseDouble(Line.split("\t")[1]) };

                //the file is later checked for correctness in when writing the values back.
                BetaBinomOverdispInSample fillerForDispersion = new BetaBinomOverdispInSample(sampleName, dispersion);
                dispersionParameters.add(fillerForDispersion);
                
            }
        }


        String binomialOutput  = FilenameUtils.getFullPath(outputLocation) + 
                                  FilenameUtils.getBaseName(outputLocation) +
                                  "_BinomialResults.txt";
        
        String betaBinomialOutput = FilenameUtils.getFullPath(outputLocation) + 
                                    FilenameUtils.getBaseName(outputLocation) +
                                    "_BetaBinomialResults.txt";
        
        
        //for ease of use initializing here.
        String CTSlinearRegressionOutput  = FilenameUtils.getFullPath(outputLocation) + 
                                  FilenameUtils.getBaseName(outputLocation) +
                                  "_CTSLinearRegressionResults.txt";
        
        
        String CTSbinomialOutput  = FilenameUtils.getFullPath(outputLocation) + 
                                  FilenameUtils.getBaseName(outputLocation) +
                                  "_CTSBinomialResults.txt";
        
        String CTSbetaBinomialOutput = FilenameUtils.getFullPath(outputLocation) + 
                                    FilenameUtils.getBaseName(outputLocation) +
                                    "_CTSBetaBinomialResults.txt";
        
        
        
        PrintWriter writer;
        int i = 0; 
        double[] dispersionVals = new double[dispersionParameters.size()]; 
        if(dispersionLocation == null){
            // use this to save the dispersionvalues
            String dispersionOutput = FilenameUtils.getFullPath(outputLocation) + 
                          FilenameUtils.getBaseName(outputLocation) +
                          "_dispersionFile.txt";
            writer = new PrintWriter(dispersionOutput, "UTF-8");       

            //header for dispersion
            writer.write("Filename\tdispersion\n");

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
        }else{

            for(BetaBinomOverdispInSample sampleDispersion : dispersionParameters){
                dispersionVals[i] = sampleDispersion.getOverdispersion()[0];
                
               
                    
                
                //do a check to make sure ordering is correct. especially important when adding a dispersion file.
                if(!(sampleDispersion.getSampleName().equals(allFiles.get(i)))){
                    System.out.println("dispersion sample name");
                    System.out.println(sampleDispersion.getSampleName());
                    System.out.println("AS file name.");
                    System.out.println(allFiles.get(i));
                    
                    throw new IllegalDataException("ERROR! ordering is not correct for filenames for overdispersion");
                }
                i++;
            }
            if(GlobalVariables.verbosity >= 10){
                System.out.println("-------------------------");
                System.out.println("Using dispersion data from a previously entered file");
                System.out.println("-------------------------");
            }
        }



        
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
                if(cellProp[i] > 1 || cellProp[i] < 0){
                    throw new IllegalDataException("The phenotype file contains values that are parsed as higher than one (1) or lower than zero (0)");
                }
                i++;
            }
        }
        
        //PART 4. Read the as files one line at a time.
        
        //Will create a file and write a header
        PrintWriter binomWriter = new PrintWriter(binomialOutput, "UTF8");
        binomWriter.println(BinomialTest.writeHeader());
        
        PrintWriter betaBinomWriter = new PrintWriter(betaBinomialOutput, "UTF8");
        betaBinomWriter.println(BetaBinomialTest.writeHeader());
        
        
        //CTS files if available.
        //CTS binomial was prone to FP, and a performance sink, therefore commented here.
        //PrintWriter CTSBinomWriter = null; 
        PrintWriter CTSBetaBinomWriter = null;
        PrintWriter CTSlinearRegressionWriter = null;
        if(phenoTypeLocation != null){
            CTSlinearRegressionWriter = new PrintWriter(CTSlinearRegressionOutput, "UTF-8");
            CTSlinearRegressionWriter.println(CTSlinearRegression.writeHeader());
            
//            CTSBinomWriter = new PrintWriter(CTSbinomialOutput, "UTF8");
//            CTSBinomWriter.println(CTSbinomialTest.writeHeader());
  
            CTSBetaBinomWriter = new PrintWriter(CTSbetaBinomialOutput, "UTF-8");
            CTSBetaBinomWriter.println(CTSBetaBinomialTest.writeHeader());
        }
        
        
        //open all the files we want to open.
        ReadAsLinesIntoIndividualSNPdata asReader = new ReadAsLinesIntoIndividualSNPdata(asLocations);

        while (true) {

            //read some the next line from the files.
            ArrayList<IndividualSnpData> allSnpData;
            allSnpData = asReader.getIndividualsFromNextLine();
            
            
            //BREAKPOINT OF THE LOOP.
            if(allSnpData.isEmpty()) break;
            
            
            
            // Add the dispersion data assuming the same ordering
            // Which was checked previously.
            for(int j = 0; j < dispersionVals.length; j++){
                allSnpData.get(j).setDispersion(dispersionVals[j]);
            }
            
            
            // add the cellProp to the snp 
            if(phenoTypeLocation != null){
                for(int j = 0; j < cellProp.length; j++){
                    allSnpData.get(j).setCellTypeProp(cellProp[j]);
                }
            }
            
            ArrayList<IndividualSnpData> het_individuals;
            het_individuals = UtilityMethods.isolateValidHeterozygotesFromIndividualSnpData(allSnpData);
    
            int numberOfHets = het_individuals.size();
            int totalOverlap = 0;

            //Data to determine If we're going to test:

            ArrayList<Integer> asRef = new ArrayList<Integer>();
            ArrayList<Integer> asAlt = new ArrayList<Integer>();
            ArrayList<Double> HetDisp = new ArrayList<Double>();
            ArrayList<Double> HetCellProp = new ArrayList<Double>();

        
            for (IndividualSnpData temp_het : het_individuals) {
                //Do nothing if there is no data in het_individuals

                asRef.add(temp_het.getRefNum());
                asAlt.add(temp_het.getAltNum());
                
                HetDisp.add(temp_het.getDispersion());
                
                if(phenoTypeLocation != null){
                    HetCellProp.add(temp_het.getCellTypeProp());
                }
                
                //this is used to check if we will continue with calculations.
                totalOverlap += temp_het.getRefNum() + temp_het.getAltNum();

            }
            
            //Print the header for the tests.
            
            if( (totalOverlap >= GlobalVariables.minReads) && (numberOfHets >= GlobalVariables.minHets) ){
            
                if(GlobalVariables.verbosity >= 10){
                    System.out.println("\n--- STARTING AS TESTS FOR: ---");
                    System.out.println("SNP name:      " + allSnpData.get(0).snpName);
                    System.out.println("at: position   " + allSnpData.get(0).chromosome + ":" + allSnpData.get(0).position + "\n");
                    System.out.println("Num of hets:      " + Integer.toString(numberOfHets));
                    
                    StringBuilder a = new StringBuilder();
                    for(Integer num : asRef){
                        a.append( String.format("% 7d ", num) );
                    }
                    
                    System.out.println("asRef:      " +  a.toString() );
                    
                    a = new StringBuilder();
                    for(Integer num : asAlt){
                        a.append( String.format("% 7d ", num) );
                    }
                    System.out.println("asAlt:      " +  a.toString() );
                    
                    a = new StringBuilder();
                    for(Double num : HetDisp){
                        a.append( String.format("% 5.4f ", num) );
                    }
                    
                    System.out.println("dispersion: " +   a.toString() );
                    
                    if(phenoTypeLocation != null){
                         a = new StringBuilder();
                        for(Double num : HetCellProp){
                            a.append( String.format("% 5.4f ", num) );
                        }
                        
                        System.out.println("cellProp:   " + a.toString() );
                    }
                }
                
                
                BinomialTest binomialResults = new BinomialTest(allSnpData);
                BetaBinomialTest betaBinomialResults = new BetaBinomialTest(allSnpData);

                if(binomialResults.isTestPerformed()){

                    binomWriter.println(binomialResults.writeTestStatistics(false));
                    betaBinomWriter.println(betaBinomialResults.writeTestStatistics(false));

                    GlobalVariables.numberOfTestPerformed++;

                }
                
                
                // do the CTS tests if data is available.

                if(phenoTypeLocation != null){
                    //do the CTS beta binomial test:
                    
                    CTSlinearRegression CTSlinearRegressionResults = new CTSlinearRegression(allSnpData);
                    
                    //the binomial test was removed because it was really slow, took about 70% of processing time, and was very false positive prone
//                    CTSbinomialTest CTSbinomialResults = new CTSbinomialTest(allSnpData, CTSlinearRegressionResults);
                    CTSBetaBinomialTest CTSbetaBinomResults = new CTSBetaBinomialTest(allSnpData, CTSlinearRegressionResults);
                    

                    // Write the results to the out_file, assuming both of them were done.
                    if (CTSbetaBinomResults.isTestPerformed()) {
                        
                       CTSlinearRegressionWriter.println(CTSlinearRegressionResults.writeTestStatistics(true));
//                       CTSBinomWriter.println(CTSbinomialResults.writeTestStatistics(true));
                       CTSBetaBinomWriter.println(CTSbetaBinomResults.writeTestStatistics(true));
                         
                    }
                }
                
                
                System.out.println("\n---- Finished SNP " + allSnpData.get(0).snpName);
            }
        }
     
        binomWriter.close();
        betaBinomWriter.close();
        
        
        if(phenoTypeLocation != null){    
            //CTSBinomWriter.close();
            CTSBetaBinomWriter.close();
        }
        
        UtilityMethods.printFinalTestStats();
        
    }


}
