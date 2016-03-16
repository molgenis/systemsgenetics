/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

import java.util.ArrayList;


/**
 *
 * @author adriaan
 */
public class BinomialTest {

    //SNP information
    private String snpName;
    private String chromosome;
    private String position;
    private String genotype;
    
    //sometimes multiple test snps have the same results.    
    private boolean TestUsedInPhasing = false;
    private ArrayList<String> additionalPositions = new ArrayList<String>();
    private ArrayList<String> additionalNames = new ArrayList<String>(); 
    
    String RegionName;
    int startOfRegion = -1;
    int endOfRegion = -1;
    int totalTestSNPs = -1;
    
    public int testRegionStart = -1;
    public int testRegionEnd = -1;
     
    
    
    
    //Information about the input data
    private int numberOfHets;
    
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
    double binomRatio;
    
    
    private boolean outPutAllData = false;
    
    
    
    
    
    // This is the constructor, will receieve an arraylist of IndividualSnpData, 
    // will filter this for only heterozygotes and will proceed if there are reads.
    // Then it will detemine the binomial proportion (no likelihood procedure, just division of ref/total).
    // Calculate log likelihoods, and do a likelihood ratio test.
    
    public BinomialTest(ArrayList<IndividualSnpData> all_individuals){
        
        //basic information, get the zero instance.
        snpName = all_individuals.get(0).getSnpName();
        chromosome = all_individuals.get(0).getChromosome();
        position = all_individuals.get(0).getPosition();
        genotype = all_individuals.get(0).genotype;
        //Isolate heterozygotes:
        ArrayList<IndividualSnpData> het_individuals = UtilityMethods.isolateValidHeterozygotesFromIndividualSnpData(all_individuals);
        
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
        
        
        if((total_overlap >= GlobalVariables.minReads) && (numberOfHets >= GlobalVariables.minHets) ){

            //Removed this because header is now done in the nonPhasedEntry and phasedEntry (to be done).

//            if(GlobalVariables.verbosity >= 10){
//
//                System.out.println("debug:");
//                System.out.println("Number of hets: " + Integer.toString(numberOfHets));
//                System.out.println(het_individuals.get(0).getSnpName());
//                System.out.println(total_overlap);
//                
//                System.out.println("asRef:          " +  asRef.toString());
//                System.out.println("asAlt:          " +  asAlt.toString());
//            }
            
            // There is data to perform the binomial test, perform it.       
            testStatistics = new BinomTest(asRef, asAlt);
            
            if(GlobalVariables.verbosity >= 10){
                System.out.println("\n--- Binomial Test Statistics: ---");

                System.out.println("\tNull logLik:   " + Double.toString(testStatistics.getNullLogLik()));
                System.out.println("\tAlt logLik:    " + Double.toString(testStatistics.getAltLogLik()));
                System.out.println();
                System.out.println("\tObserved Prop: " + Double.toString(testStatistics.getBinomRatio()));
                System.out.println();
                System.out.println("\tchi-sq:        " + Double.toString(testStatistics.getChiSq()));
                System.out.println("\tP value:       " + Double.toString(testStatistics.getpVal()));
                System.out.println("--------------------------------");
            }
            
            //binomial_test perfomed, will now set test performed to true
            testPerformed = true;
            binomRatio = testStatistics.getBinomRatio();
        } else{
            // there is no data to do the binomial test,
            // Will not set variables.
        }
    }
    
    //constructor method didn't work, so doing it like this.
    public static  BinomialTest phasedBinomialTest(ArrayList<IndividualSnpData> all_individuals, GenomicRegion thisRegion, int testSNPs){
        
        BinomialTest t = new BinomialTest(all_individuals);
        
        
        t.TestUsedInPhasing = true;
        t.RegionName = thisRegion.getAnnotation();
        t.startOfRegion = thisRegion.getStartPosition();
        t.endOfRegion = thisRegion.getEndPosition();
        
        t.testRegionStart = thisRegion.getTestStart();
        t.testRegionEnd   = thisRegion.getTestEnd();
        
        
        t.totalTestSNPs = testSNPs;
        
        return t;
    
    } 
    
   public static String writeHeader(){
       String header = "chr\tpos\tsnpName\tnumHets\tpVal\tchiSq\tbinomRatio\tnullLogLik\taltLogLik";
       return header;
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

            


        } else {
            //when no testing is done, will only output snp name and position, and NA.
            for(int i=0; i < 5; i++ ){
                out += "NA\t";
            
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
    

    /**
     * @return the multipleTestSnps
     */
    public boolean hasMultipleTestSnps() {
        return TestUsedInPhasing;
    }

    /**
     * @return the additionalPositions
     */
    public ArrayList<String> getAdditionalPositions() {
        return additionalPositions;
    }

    /**
     * @return the additionalNames
     */
    public ArrayList<String> getAdditionalNames() {
        return additionalNames;
    }

    public void addAdditionalSNP(String snpName, String snpPos){
        additionalNames.add(snpName);
        additionalPositions.add(snpPos);
    
    }

    void setSnpName(String snpName1) {
        snpName = snpName1 ;
    }

    /**
     * @return the genotype
     */
    public String getGenotype() {
        return genotype;
    }

    /**
     * @param genotype the genotype to set
     */
    public void setGenotype(String genotype) {
        this.genotype = genotype;
    }
    
}

