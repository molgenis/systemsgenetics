/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cell_type_specific_ase;

import java.util.ArrayList;


/**
 *
 * @author adriaan
 */
public class BinomialTest {

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
    
    //Output of the actual test.
    
    private BinomTest testStatistics;
    private boolean testPerformed = false; 
    
    
    private final boolean outPutAllData = false;
    
    // This is the constructor, will receieve an arraylist of IndividualSnpData, 
    // will filter this for only heterozygotes and will proceed if there are reads.
    // Then it will detemine the binomial proportion (no likelihood procedure, just division of ref/total).
    // Calculate log likelihoods, and do a likelihood ratio test.
    
    public BinomialTest(ArrayList<IndividualSnpData> all_individuals){
        boolean debug = true;
        //basic information, get the zero instance.
        snpName = all_individuals.get(0).getSnpName();
        chromosome = all_individuals.get(0).getChromosome();
        position = all_individuals.get(0).getPosition();
        
        //Isolate heterozygotes:
        ArrayList<IndividualSnpData> het_individuals = all_individuals;
        
        numberOfHets = het_individuals.size();
        
        hetSampleNames = new ArrayList<String>();
        asRef = new ArrayList<Integer>(); 
        asAlt = new ArrayList<Integer>(); 
        asNo  = new ArrayList<Integer>(); 
        
        
        int total_overlap = 0;
        
        for (IndividualSnpData temp_het : het_individuals) {
            //Do nothing if there is no data in het_individuals

            hetSampleNames.add(temp_het.getSampleName());
            
            asRef.add(temp_het.getRefNum());
            asAlt.add(temp_het.getAltNum());
            asNo.add(temp_het.getNoNum());
            
            //this is used to check if we will continue with calculations.
            total_overlap += temp_het.getRefNum() + temp_het.getAltNum();
        }
        
        
        if(total_overlap > 10){
            if(debug){
                System.out.println("debug:");
                System.out.println("Number of hets: " + Integer.toString(numberOfHets));
                System.out.println(het_individuals.get(0).getSnpName());
                System.out.println(total_overlap);
                
                System.out.println("asRef:" +  asRef.toString());
                System.out.println("asAlt:" +  asAlt.toString());
            }
            // There is data to perform the binomial test, perform it.       
            testStatistics = new BinomTest(asRef, asAlt);
            //binomial_test perfomed, will now set test performed to true
            testPerformed = true;
            
        } else{
            // there is no data to do the binomial test,
            // Will not set variables.
            // Don't know what to do here, should I initialize a BinomTestStatistic Object?.    
        }
    }
    
    private ArrayList<IndividualSnpData> isolateHeterozygotes(ArrayList<IndividualSnpData> all_individuals) {
        
        ArrayList<IndividualSnpData> hets;
        hets = new ArrayList<IndividualSnpData>();
        
        for(IndividualSnpData sample : all_individuals){
            
            String genotype = sample.getGenotype();
            
            //assuming the genotype is formatted as: "[C, A]"
            
            char charOne = genotype.charAt(1);
            char charTwo = genotype.charAt(4);
            
            if(charOne != charTwo){
                hets.add(sample);
            }       
        }
        
        return hets;
    }

   public String writeTestStatistics(boolean outPutAllData){
        
        String out = "";
        
        //snp info
        out += chromosome + "\t";
        out += position + "\t";
        out += snpName + "\t";
        
        if(testPerformed){
            out += numberOfHets + "\t";

            // lets have a look at how the decimals behave here,
            // otherwise i think I have to do something with printf or something.
            out += Double.toString(testStatistics.getpVal()) + "\t";
            out += Double.toString(testStatistics.getChiSq()) + "\t";
            out += Double.toString(testStatistics.getBinomRatio()) + "\t";
            out += Double.toString(testStatistics.getNullLogLik()) + "\t";
            out += Double.toString(testStatistics.getAltLogLik());

            if(outPutAllData){
                
                String samples_string="";
                String ref_string="";
                String alt_string="";
                String no_string="";

                for(int i=0; i < hetSampleNames.size(); i++){
                    
                    //samples_string += hetSampleNames.get(i) + ";";
                    ref_string += Integer.toString(asRef.get(i)) + ";";
                    alt_string += Integer.toString(asAlt.get(i)) + ";";
                    no_string += Integer.toString(asNo.get(i)) + ";";
                }
                
                //remove last delimiter
                //samples_string = samples_string.substring(0, samples_string.length()-1);
                ref_string = ref_string.substring(0, ref_string.length()-1);
                alt_string = alt_string.substring(0, alt_string.length()-1);
                no_string = no_string.substring(0, no_string.length()-1);
                
                //out += "\t" + samples_string + "\t" + ref_string + "\t" + alt_string + "\t" + no_string;
                out += "\t" + ref_string + "\t" + alt_string + "\t" + no_string;


            }


        } else {
            //when no testing is done, will only output snp name and position, and NA.
            for(int i=0; i < 5; i++ ){
                out += "NA\t";
            
            }
            if(outPutAllData){
                for(int i=0; i < 4; i++ ){
                    out += "NA\t";

                }
            }
            out += "NA";
        
        }
        
        return out;
    }
    
    
    /**
     * @return the snpName
     */
    public String getSnpName() {
        return snpName;
    }

    /**
     * @return the chromosome
     */
    public String getChromosome() {
        return chromosome;
    }

    /**
     * @return the position
     */
    public String getPosition() {
        return position;
    }

    /**
     * @return the numberOfHets
     */
    public int getNumberOfHets() {
        return numberOfHets;
    }

    /**
     * @return the hetSampleNames
     */
    public ArrayList<String> getHetSampleNames() {
        return hetSampleNames;
    }

    /**
     * @return the asRef
     */
    public ArrayList<Integer> getAsRef() {
        return asRef;
    }

    /**
     * @return the asAlt
     */
    public ArrayList<Integer> getAsAlt() {
        return asAlt;
    }

    /**
     * @return the asNo
     */
    public ArrayList<Integer> getAsNo() {
        return asNo;
    }

    /**
     * @return the testStatistics
     */
    public BinomTest getTestStatistics() {
        return testStatistics;
    }

    /**
     * @return the testPerformed
     */
    public boolean isTestPerformed() {
        return testPerformed;
    }
    
    
    
}

