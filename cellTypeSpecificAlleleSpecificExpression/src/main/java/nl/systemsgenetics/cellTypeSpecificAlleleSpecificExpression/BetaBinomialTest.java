/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

import java.util.ArrayList;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;

/**
 *
 * @author adriaan
 */

class BetaBinomialTest {
    
    //SNP information
    private final String snpName;
    private final String chromosome;
    private final String position;
    private String genotype;
    
    //Information about the input data
    private final int numberOfHets;
    
    //sometimes multiple test snps have the same results, this is done 
    //deduplication purposes.
    
    boolean TestUsedInPhasing = false;
    private ArrayList<String> additionalPositions = new ArrayList<String>();
    private ArrayList<String> additionalNames = new ArrayList<String>(); 
    
    String RegionName;
    int startOfRegion = -1;
    int endOfRegion = -1;
        
    public int testRegionStart = -1;
    public int testRegionEnd = -1;
    
    
    int totalTestSNPs = -1;
    
    //Names of individuals for easy reference
    
    //The ArrayLists below are for 
    private ArrayList<String> hetSampleNames;
    
    //Data that will be actually used
    
    private ArrayList<Integer> asRef;
    private ArrayList<Integer> asAlt;
    private ArrayList<Integer> asNo;
    
    private ArrayList<Double> dispersion;
    
    //Output of the actual test.
    private boolean testPerformed = false;
    private final boolean outPutAllData = false;
    
    //Test statistics:
    
    //MLE statistics
    int iterations;
    double alphaParam;
    double betaParam;
    double binomRatio;
    
    //likelihood of data:
    double nullLogLik;
    double altLogLik;
    double chiSq;
    double pVal;
    

    
    public BetaBinomialTest(ArrayList<IndividualSnpData> all_individuals) throws Exception{
        boolean debug=true;
        //basic information, get the zero instance, was checked as the same in ReadAsLinesIntoIndividualSNPdata
        snpName = all_individuals.get(0).getSnpName();
        chromosome = all_individuals.get(0).getChromosome();
        position = all_individuals.get(0).getPosition();
        genotype = all_individuals.get(0).genotype;
        
        ArrayList<IndividualSnpData> het_individuals;
        het_individuals = UtilityMethods.isolateValidHeterozygotesFromIndividualSnpData(all_individuals);
    
        numberOfHets = het_individuals.size();
        
        hetSampleNames = new ArrayList<String>();
        asRef = new ArrayList<Integer>(); 
        asAlt = new ArrayList<Integer>(); 
        asNo  = new ArrayList<Integer>(); 
        
        dispersion =  new ArrayList<Double>(); 
        
        int total_overlap = 0;
        

        
        
        for (IndividualSnpData temp_het : het_individuals) {
            //Do nothing if there is no data in het_individuals

            hetSampleNames.add(temp_het.getSampleName());
            
            asRef.add(temp_het.getRefNum());
            asAlt.add(temp_het.getAltNum());
            asNo.add(temp_het.getNoNum());
            
            dispersion.add(temp_het.getDispersion());
            
            
            //this is used to check if we will continue with calculations.
            total_overlap += temp_het.getRefNum() + temp_het.getAltNum();
            
        }
        
        
        if((total_overlap >= GlobalVariables.minReads) && (numberOfHets >= GlobalVariables.minHets) ){
            // There is data to perform the binomial test, perform it.       
            if(GlobalVariables.verbosity >= 100){
                System.out.println();
                System.out.println("---- Starting beta binomial LRT test estimate ----");
                System.out.println("SNP name: " + snpName);
                System.out.println("at: chr" + chromosome + ":" + position);
                System.out.println("--------------------------------------------------");

                
                System.out.println("debug:");
                System.out.println("Num of hets: " + Integer.toString(numberOfHets));
                System.out.println(het_individuals.get(0).getSnpName());
                System.out.println(total_overlap);
                
                System.out.println("asRef:       " +  asRef.toString());
                System.out.println("asAlt:       " +  asAlt.toString());
                System.out.println("dispersion:  " +  dispersion.toString());
                System.out.println("Starting Null estimation");
            }

            
            
            Integer[] asRefArray = asRef.toArray(new Integer[asRef.size()]);
            Integer[] asAltArray = asAlt.toArray(new Integer[asAlt.size()]);
            Double[]  dispArray  = dispersion.toArray(new Double[dispersion.size()]);
            
            
            
            BetaBinomNullLikelihood betaBinomNull;
            betaBinomNull = new BetaBinomNullLikelihood(asRefArray, 
                                                        asAltArray, 
                                                        dispArray
                                                        );
            
            
            nullLogLik = betaBinomNull.value(new double[] {0.5});
            

            if(GlobalVariables.verbosity >= 100){
                System.out.println("Starting Alt estimation");
            }
            
            BetaBinomAltLikelihood betaBinomAlt;
            betaBinomAlt = new BetaBinomAltLikelihood(asRefArray, 
                                                        asAltArray, 
                                                        dispArray
                                                        );
            NelderMeadSimplex simplex;
            simplex = new NelderMeadSimplex(2);
            SimplexOptimizer optimizer = new SimplexOptimizer(GlobalVariables.simplexThreshold, GlobalVariables.simplexThreshold); //numbers are to which precision you want it to be done.
            PointValuePair solutionAlt = optimizer.optimize(
                                            new ObjectiveFunction(betaBinomAlt),
                                            new MaxEval(GlobalVariables.maximumIterations),
                                            simplex,
                                            GoalType.MINIMIZE,
                                            new InitialGuess(new double[] {0.5, 0.5}), 
                                            new SearchInterval(-1000.0, 1000.0)
                                            );
            
            double[] valueAlt = solutionAlt.getPoint();
            
            altLogLik  = betaBinomAlt.value(valueAlt);
            iterations = optimizer.getIterations();
            alphaParam = valueAlt[0];
            betaParam = valueAlt[1];
            binomRatio = valueAlt[0] / (valueAlt[0] + valueAlt[1]);
            
            //chi squared statistic is determined based on both null and alt loglikelihoods.
            chiSq = 2.0 * (nullLogLik - altLogLik);

            //determine P value based on distribution
            ChiSquaredDistribution distribution = new ChiSquaredDistribution(1);
            pVal = 1 - distribution.cumulativeProbability(chiSq);
            
            if(GlobalVariables.verbosity >= 10){
                System.out.println("\n--- Beta Binomial Test Statistics: ---");
                System.out.println("LogLik of Alt converged to a threshold of " + Double.toString(GlobalVariables.simplexThreshold));
                System.out.println("\tAlpha parameter:      " + Double.toString(valueAlt[0]));
                System.out.println("\tBeta parameter:       " + Double.toString(valueAlt[1]));
                System.out.println("\tIterations to converge:           " + Integer.toString(iterations) + "\n");
                System.out.println("\tNull log likelihood:  " + Double.toString(nullLogLik));   
                System.out.println("\tAlt log likelihood:   " + Double.toString(altLogLik) + "\n");
                System.out.println("\tChisq statistic:      " + Double.toString(chiSq));
                System.out.println("\tP value:              " + Double.toString(pVal));
                System.out.println("------------------------------------------");
                //TODO, I want to format this properly, but not necessary

            }
            
            testPerformed = true;
            
            //This below is to compare everything to python created by WASP
            //Will put the backup back into python, just to make sure.
            
            if(GlobalVariables.verbosity >= 100){
                System.out.println("loglikRef = " + Double.toString(altLogLik));
                System.out.println("loglik = 0");
                for(int i = 0; i < asRef.size();i++){

                    System.out.println("loglik += AS_betabinom_loglike("
                                        + "[math.log(" + Double.toString(alphaParam) +") - math.log(" + Double.toString(alphaParam) + " + " + Double.toString(betaParam) +"), " +
                                       "math.log(" + Double.toString(betaParam) +") - math.log(" + Double.toString(alphaParam) + " + " + Double.toString(betaParam) + ")],"
                                       + Double.toString(dispersion.get(i)) + ","
                                       + Double.toString(asRef.get(i)) + ","
                                       + Double.toString(asAlt.get(i)) + ","   
                                       + "0.980198, 0.005)");
                }
                System.out.println("if allclose(-1.0 * loglik, loglikRef): correctTest += 1");
            

                //I find exactly the correct values after running it through WASP implementation in python, after one check.
                //Using this as a foundation stone to write some tests:
                for(int i = 0; i < asRef.size();i++){
                    String command = "Assert.assertEquals( " +
                             "likelihoodFunctions.BetaBinomLogLik(" + String.format("%.15g", dispArray[i]) + ","
                              + String.format("%.15g", alphaParam) + ","
                              + String.format("%.15g", betaParam) + "," 
                              + Integer.toString( asRefArray[i]) + ","
                             + Integer.toString(asAltArray[i]) + ")" +
                            ", " + Double.toString(LikelihoodFunctions.BetaBinomLogLik(dispArray[i], alphaParam, betaParam, asRefArray[i], asAltArray[i])) + ");" ;
                     System.out.println(command);
                }
            }
        
        }
    }

    //constructor method didn't work, so doing it like this.
    public static  BetaBinomialTest phasedBetaBinomialTest(ArrayList<IndividualSnpData> all_individuals, GenomicRegion thisRegion, int testSNPs) throws Exception{
        
        BetaBinomialTest t = new BetaBinomialTest(all_individuals);
        
        
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
       String header = "chr\tpos\tsnpName\tnumHets\tpVal\tchiSq\tbinomRatio\tnullLogLik\taltLogLik\talphaParam\tbetaParam";
       return header;
    } 
    

    public String writeTestStatistics(boolean all_data) {
        
        String out = "";
        
        //snp info
        out += chromosome + "\t";
        out += position + "\t";
        out += snpName + "\t";
        
        if(testPerformed){
            out += numberOfHets + "\t";

            // lets have a look at how the decimals behave here, seems all right.
            // otherwise i think I have to do something with printf or something.
            out += Double.toString(pVal) + "\t";
            out += Double.toString(chiSq) + "\t";
            out += Double.toString(binomRatio) + "\t";
            out += Double.toString(nullLogLik) + "\t";
            out += Double.toString(altLogLik) + "\t";
            out += Double.toString(alphaParam) + "\t";
            out += Double.toString(betaParam);
            
            
        } else {
            //when no testing is done, will only output snp name and position, and NA.
            //Make sure this is still correct.
            
            for(int i=0; i < 6; i++ ){
                out += "NA\t";
            
            }

            out += "NA";
        
        }
        return out;   
    }
    
    
    
    /**
     * @return the testPerformed
     */
    public boolean isTestPerformed() {
        return testPerformed;
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
