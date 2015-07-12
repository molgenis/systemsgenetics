/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

import java.util.ArrayList;
import static java.util.Collections.list;
import java.util.List;
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

    //Information about the input data
    private final int numberOfHets;
    
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
    
    
    private final double precision = 1.0E-7;  
    
    public BetaBinomialTest(ArrayList<dispersedIndividualSnpData> all_individuals, int minReads, int minHets){
        boolean debug=true;
        //basic information, get the zero instance.
        snpName = all_individuals.get(0).getSnpName();
        chromosome = all_individuals.get(0).getChromosome();
        position = all_individuals.get(0).getPosition();
        
        //isolate heterozygote individuals.
        ArrayList<dispersedIndividualSnpData> het_individuals = isolateHeterozygotes(all_individuals);
    
        numberOfHets = het_individuals.size();
        
        hetSampleNames = new ArrayList<String>();
        asRef = new ArrayList<Integer>(); 
        asAlt = new ArrayList<Integer>(); 
        asNo  = new ArrayList<Integer>(); 
        
        dispersion =  new ArrayList<Double>(); 
        
        int total_overlap = 0;
        

        
        
        for (dispersedIndividualSnpData temp_het : het_individuals) {
            //Do nothing if there is no data in het_individuals

            hetSampleNames.add(temp_het.getSampleName());
            
            asRef.add(temp_het.getRefNum());
            asAlt.add(temp_het.getAltNum());
            asNo.add(temp_het.getNoNum());
            
            dispersion.add(temp_het.getDispersion());
            
            
            //this is used to check if we will continue with calculations.
            total_overlap += temp_het.getRefNum() + temp_het.getAltNum();
            
        }
        
        
        if((total_overlap >= minReads) && (numberOfHets >= minHets) ){
            // There is data to perform the binomial test, perform it.       
            System.out.println();
            System.out.println("---- Starting beta binomial LRT test estimate ----");
            System.out.println("SNP name: " + snpName);
            System.out.println("at: chr" + chromosome + ":" + position);
            System.out.println("--------------------------------------------------");
            
            if(debug){
                
                System.out.println("debug:");
                System.out.println("Num of hets: " + Integer.toString(numberOfHets));
                System.out.println(het_individuals.get(0).getSnpName());
                System.out.println(total_overlap);
                
                System.out.println("asRef:       " +  asRef.toString());
                System.out.println("asAlt:       " +  asAlt.toString());
                System.out.println("dispersion:  " +  dispersion.toString());
            
            }

            
            
            Integer[] asRefArray = asRef.toArray(new Integer[asRef.size()]);
            Integer[] asAltArray = asAlt.toArray(new Integer[asAlt.size()]);
            Double[]  dispArray  = dispersion.toArray(new Double[dispersion.size()]);
            
            System.out.println("Starting Null estimation");
            
            BetaBinomNullLikelihood betaBinomNull;
            betaBinomNull = new BetaBinomNullLikelihood(asRefArray, 
                                                        asAltArray, 
                                                        dispArray
                                                        );
            
            
            nullLogLik = betaBinomNull.value(new double[] {0.5});
            

            
            System.out.println("Starting Alt estimation");

            
            BetaBinomAltLikelihood betaBinomAlt;
            betaBinomAlt = new BetaBinomAltLikelihood(asRefArray, 
                                                        asAltArray, 
                                                        dispArray
                                                        );
            NelderMeadSimplex simplex;
            simplex = new NelderMeadSimplex(2);
            SimplexOptimizer optimizer = new SimplexOptimizer(precision, precision); //numbers are to which precision you want it to be done.
            PointValuePair solutionAlt = optimizer.optimize(
                                            new ObjectiveFunction(betaBinomAlt),
                                            new MaxEval(500),
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
            
            
            System.out.println("LogLik of Alt converged to a threshold of " + Double.toString(precision));
            System.out.println("\tAlpha parameter:      " + Double.toString(valueAlt[0]));
            System.out.println("\tBeta parameter:       " + Double.toString(valueAlt[1]));
            System.out.println("\tIterations to converge:           " + Integer.toString(iterations) + "\n");
            System.out.println("\tNull log likelihood:  " + Double.toString(nullLogLik));   
            System.out.println("\tAlt log likelihood:   " + Double.toString(altLogLik) + "\n");
            System.out.println("\tChisq statistic:      " + Double.toString(chiSq));
            System.out.println("\tP value:              " + Double.toString(pVal));
            //TODO, I want to format this properly, but not necessary
            System.out.println("\n---- Finished SNP " + snpName);
            testPerformed = true;
            
        }
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
            
            if(outPutAllData){
                String samples_string="";
                String ref_string="";
                String alt_string="";                
                String no_string="";
                String disp_string="";
                
                for(int i=0; i < hetSampleNames.size(); i++){
                    
                    //samples_string += hetSampleNames.get(i) + ";";
                    ref_string += Integer.toString(asRef.get(i)) + ";";
                    alt_string += Integer.toString(asAlt.get(i)) + ";";
                    no_string += Integer.toString(asNo.get(i)) + ";";
                    disp_string += Double.toString(dispersion.get(i)) + ";";
                }
                
                //remove last delimiter
                //samples_string = samples_string.substring(0, samples_string.length()-1);
                ref_string = ref_string.substring(0, ref_string.length()-1);
                alt_string = alt_string.substring(0, alt_string.length()-1);
                no_string = no_string.substring(0, no_string.length()-1);
                disp_string = disp_string.substring(0, disp_string.length()-1);
                
                //out += "\t" + samples_string + "\t" + ref_string + "\t" + alt_string + "\t" + no_string;
                out += "\t" + ref_string + "\t" + alt_string + "\t" + no_string  + "\t" + disp_string;


            }


        } else {
            //when no testing is done, will only output snp name and position, and NA.
            //Make sure this is still correct.
            
            for(int i=0; i < 6; i++ ){
                out += "NA\t";
            
            }
            if(outPutAllData){
                for(int i=0; i < 5; i++ ){
                    out += "NA\t";

                }
            }
            out += "NA";
        
        }
        return out;   
    }

    
    private ArrayList<dispersedIndividualSnpData> isolateHeterozygotes(ArrayList<dispersedIndividualSnpData> all_individuals) {
        
        ArrayList<dispersedIndividualSnpData> hets;
        hets = new ArrayList<dispersedIndividualSnpData>();
        
        for(dispersedIndividualSnpData sample : all_individuals){
            
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

}