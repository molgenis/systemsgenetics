/*
 * Copyright (C) 2015 adriaan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

import java.util.ArrayList;
import org.apache.commons.math3.stat.regression.RegressionResults;
import org.apache.commons.math3.stat.regression.SimpleRegression;

/**
 *
 * @author adriaan
 */
class CTSlinearRegression {

    //SNP information
    private final String snpName;
    private final String chromosome;
    private final String position;

    //Information about the input data
    private final int numberOfHets;
    
    //Names of individuals for easy reference
    
    //The ArrayLists below are for 
    private ArrayList<String> hetSampleNames;
    
    //Data that will be actually used
    
    private ArrayList<Integer> asRef;
    private ArrayList<Integer> asAlt;
    private ArrayList<Integer> asNo;
    
    private ArrayList<Double> cellProp;
    
    boolean testPerformed = false;
    private final boolean outPutAllData = false;
    
    Double pValue;
    Double slope = 0.0;
    Double intercept = 0.0;
    
    
    Double stdErrorSlope;
    Double stdErrorIntercept;
    
    Double Rsquared;
    
    
    public CTSlinearRegression(ArrayList<IndividualSnpData> all_individuals) {
       
        //basic information, get the zero instance.
        snpName = all_individuals.get(0).getSnpName();
        chromosome = all_individuals.get(0).getChromosome();
        position = all_individuals.get(0).getPosition();
        
        //isolate heterozygotes
        ArrayList<IndividualSnpData> het_individuals = UtilityMethods.isolateValidHeterozygotesFromIndividualSnpData(all_individuals);
        numberOfHets = het_individuals.size();
        
                hetSampleNames = new ArrayList<String>();
        asRef = new ArrayList<Integer>(); 
        asAlt = new ArrayList<Integer>(); 
        asNo  = new ArrayList<Integer>(); 
        
      
        cellProp   = new ArrayList<Double>();
        int total_overlap = 0;
        
        //Get the basic data without doing any tests.
        for (IndividualSnpData temp_het : het_individuals) {
            //Do nothing if there is no data in het_individuals

            hetSampleNames.add(temp_het.getSampleName());
            
            asRef.add(temp_het.getRefNum());
            asAlt.add(temp_het.getAltNum());
            asNo.add(temp_het.getNoNum());
            cellProp.add(temp_het.getCellTypeProp());
            
            //this is used to check if we will continue with calculations.
            //BASED on the minHets and MinReads
            total_overlap += temp_het.getRefNum() + temp_het.getAltNum();
        }
        
        
        //Check if we do a test.
        if((total_overlap >= GlobalVariables.minReads) && (numberOfHets >= GlobalVariables.minHets) && (numberOfHets >= 3)){
            
            ASScatterPlot plotThis = null;
            
            if(!GlobalVariables.plotDir.equals("")){
                plotThis = new ASScatterPlot(400);
            }
            
            
            SimpleRegression thisRegression = new SimpleRegression();
            for(int i=0; i< asRef.size(); i++ ){
                Double asRatio;
                //do this check, otherwise the denominator will be zero.
                
                if(asRef.get(i) != 0){
                    asRatio = ((double)asRef.get(i)) / ((double)(asRef.get(i) + asAlt.get(i)));
                }else{
                    asRatio = 0.0;
                }
                
                Double phenoRatio =  cellProp.get(i);
                thisRegression.addData(phenoRatio, asRatio);
                if(!GlobalVariables.plotDir.equals("")){
                    plotThis.plot(asRatio, phenoRatio);
                }
            }
            
            if(!GlobalVariables.plotDir.equals("")){
                plotThis.draw(GlobalVariables.plotDir + "/" + snpName +  "_ASratio_Pheno_Plot.png");
            } 
            
            slope = thisRegression.getSlope();
            intercept  = thisRegression.getIntercept();
            Rsquared = thisRegression.getRSquare();
            stdErrorIntercept = thisRegression.getInterceptStdErr();
            stdErrorSlope = thisRegression.getSlopeStdErr();
            
            pValue = thisRegression.getSignificance();
            
            
            if(GlobalVariables.verbosity >= 10){
                    System.out.println("\n--- Starting cell type specific linear regression ---");
                    System.out.println("\tSlope:                   " + Double.toString(slope));
                    System.out.println("\tStdError of Slope:       " + Double.toString(stdErrorSlope) + "\n");
                    System.out.println("\tIntercept:               " + Double.toString(intercept));
                    System.out.println("\tStdError of Intercept:   " + Double.toString(stdErrorIntercept) + "\n");   
                    System.out.println("\tP value:                 " + Double.toString(pValue));
                    System.out.println("--------------------------------------------------------------");
                    
            }
            
            
            testPerformed = true;
        }
        
    }
    
    public static String writeHeader(){
       String header = "chr\tpos\tsnpName\tnumHets\tpVal\tslope\tintercept\tstdErrorSlope\tstdErrorIntercept";
       return header;
    } 
    
    
    public String writeTestStatistics(boolean output_all_data){
        
        String out = "";
        
        //snp info
        out += chromosome + "\t";
        out += position + "\t";
        out += snpName + "\t";
        
        if(testPerformed){
            out += numberOfHets + "\t";
            
            out += pValue + "\t";
            out += slope + "\t";
            out += intercept + "\t";
            
            out +=  stdErrorIntercept + "\t";
            out +=  stdErrorSlope;  
            
            
        } else {
            
            for(int i=0; i < 5; i++ ){
                out += "NA\t";
            
            }
            
            out += "NA";            
        }
        
        return out;
        
    }
    
}
