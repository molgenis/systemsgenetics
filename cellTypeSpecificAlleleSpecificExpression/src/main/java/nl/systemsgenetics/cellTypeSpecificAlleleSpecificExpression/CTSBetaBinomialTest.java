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
class CTSBetaBinomialTest {

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
    private ArrayList<Double> cellProp;
    
    
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
    
    //precision
    double precision = 1E-6;
    
    public CTSBetaBinomialTest(ArrayList<CTSdispersedIndividualSnpData> all_individuals, int minReads, int minHets) {
        
        CTSdispersedIndividualSnpData
    
        boolean debug=true;
        //basic information, get the zero instance.
        snpName = all_individuals.get(0).getSnpName();
        chromosome = all_individuals.get(0).getChromosome();
        position = all_individuals.get(0).getPosition();
        
        //isolate heterozygote individuals.
        ArrayList<CTSdispersedIndividualSnpData> het_individuals = isolateHeterozygotes(all_individuals);
    
        numberOfHets = het_individuals.size();
        
        hetSampleNames = new ArrayList<String>();
        asRef = new ArrayList<Integer>(); 
        asAlt = new ArrayList<Integer>(); 
        asNo  = new ArrayList<Integer>(); 
        
        dispersion =  new ArrayList<Double>(); 
        
        int total_overlap = 0;
        
        //Get the basic data without doing any tests.
        for (CTSdispersedIndividualSnpData temp_het : het_individuals) {
            //Do nothing if there is no data in het_individuals

            hetSampleNames.add(temp_het.getSampleName());
            
            asRef.add(temp_het.getRefNum());
            asAlt.add(temp_het.getAltNum());
            asNo.add(temp_het.getNoNum());
            
            dispersion.add(temp_het.getDispersion());
            
            
            //this is used to check if we will continue with calculations.
            //BASED on the minHets and MinReads
            total_overlap += temp_het.getRefNum() + temp_het.getAltNum();
        }
        
        //Check if we do a test.
        if((total_overlap >= minReads) && (numberOfHets >= minHets) ){

            // There is data to perform the binomial test, perform it.       
            System.out.println();
            System.out.println("---- Starting CTS beta binomial LRT test estimate ----");
            System.out.println("SNP name: " + snpName);
            System.out.println("at: chr" + chromosome + ":" + position);
            System.out.println("------------------------------------------------------");
            
            if(debug){
                
                System.out.println("debug:");
                System.out.println("Num of hets: " + Integer.toString(numberOfHets));
                System.out.println(het_individuals.get(0).getSnpName());
                System.out.println(total_overlap);
                
                System.out.println("asRef:       " +  asRef.toString());
                System.out.println("asAlt:       " +  asAlt.toString());
                System.out.println("dispersion:  " +  dispersion.toString());
            
            }
            
            /*
                Going to do this in a two step procedure:
                
                1. Determine a   2 parameter null the same way as what is done 
                        in the beta binomial alt, with starting values:
                        This will be the actual null model in the CTS beta binomial output
                2. Determine the cell type specific proportion, but now we use a
                        cell type specific model as the alternative.
                        
            
            */
            
            Integer[] asRefArray      = asRef.toArray(new Integer[asRef.size()]);
            Integer[] asAltArray      = asAlt.toArray(new Integer[asAlt.size()]);
            Double[]  dispArray       = dispersion.toArray(new Double[dispersion.size()]);
            Double[]  cellPropArray   = cellProp.toArray(new Double[cellProp.size()]);
            
            
            System.out.println("Starting non-CTS Beta binomial Null estimation");
            
            
            BetaBinomAltLikelihood betaBinomNull;
            betaBinomNull = new BetaBinomAltLikelihood(asRefArray, 
                                                        asAltArray, 
                                                        dispArray
                                                        );
            NelderMeadSimplex simplex;
            simplex = new NelderMeadSimplex(2);
            SimplexOptimizer optimizer = new SimplexOptimizer(precision, precision); //numbers are to which precision you want it to be done.
            PointValuePair solutionNull = optimizer.optimize(
                                            new ObjectiveFunction(betaBinomNull),
                                            new MaxEval(500),
                                            simplex,
                                            GoalType.MINIMIZE,
                                            new InitialGuess(new double[] {0.5, 0.5}), 
                                            new SearchInterval(-1000.0, 1000.0)
                                            );
            
            double[] valueNull = solutionNull.getPoint();
            
            altLogLik  = betaBinomNull.value(valueNull);
            iterations = optimizer.getIterations();
            alphaParam = valueNull[0];
            betaParam = valueNull[1];
            binomRatio = valueNull[0] / (valueNull[0] + valueNull[1]);
            
            System.out.println("LogLik of Null converged to a threshold of " + Double.toString(precision));
            System.out.println("\tAlpha parameter:      " + Double.toString(valueNull[0]));
            System.out.println("\tBeta parameter:       " + Double.toString(valueNull[1]));
            System.out.println("\tIterations to converge:           " + Integer.toString(iterations) + "\n");
            
            
            CTSbetaBinomialAltLikelihood CTSbetaBinomAlt;
            CTSbetaBinomAlt = new CTSbetaBinomialAltLikelihood(asRefArray, 
                                                        asAltArray, 
                                                        cellPropArray,
                                                        dispArray
                                                        );
            
            simplex = new NelderMeadSimplex(4);
            PointValuePair solutionAlt = optimizer.optimize(
                                            new ObjectiveFunction(CTSbetaBinomAlt),
                                            new MaxEval(500),
                                            simplex,
                                            GoalType.MINIMIZE,
                                            new InitialGuess(new double[] {0.0, 0.0,valueNull[0], valueNull[1] }), //Start with the loglik of the null.
                                            new SearchInterval(-1000.0, 1000.0)
                                            );
            
            double[] valueAlt = solutionAlt.getPoint();
            
            
            altLogLik  = CTSbetaBinomAlt.value(valueAlt);
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
            
            
            
            
            //Finally this was done, so we say the test was performed.
            testPerformed = true;
            
        }
        
    
    }
    
    
    
    private ArrayList<CTSdispersedIndividualSnpData> isolateHeterozygotes(ArrayList<CTSdispersedIndividualSnpData> all_individuals) {
        
        ArrayList<CTSdispersedIndividualSnpData> hets;
        hets = new ArrayList<CTSdispersedIndividualSnpData>();
        
        for(CTSdispersedIndividualSnpData sample : all_individuals){
            
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
    
    
    
    
}
