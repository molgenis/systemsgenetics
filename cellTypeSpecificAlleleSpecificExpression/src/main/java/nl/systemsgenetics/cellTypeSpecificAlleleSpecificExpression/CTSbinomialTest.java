/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

import java.util.ArrayList;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
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
class CTSbinomialTest {
    
    //SNP information
    private final String snpName;
    private final String chromosome;
    private final String position;

    //Information about the input data
    private final int numberOfHets;
    
    //Names of individuals for easy reference
    
     
    private ArrayList<String> hetSampleNames;
    
    //The ArrayLists below are for
    //Data that will be actually used
    //In calculation.
    
    private ArrayList<Integer> asRef;
    private ArrayList<Integer> asAlt;
    private ArrayList<Integer> asNo;
    
    private ArrayList<Double> cellProp;
    
    
    
    
    
    //Specific parameters for the cell type specific test
    double MLEratioCellType;
    double MLEratioResidual;
    
    
    
    //standard test output:    
    int iterations;
    double nullLogLik;
    double altLogLik;
    double pVal;
    double chiSq;
    
    
    //some settings
    boolean outPutAllData = false;
    boolean testPerformed = false;
    boolean testConverged  = false;
    
    
    public CTSbinomialTest(ArrayList<IndividualSnpData> all_individuals, CTSlinearRegression CTSlinearRegressionResults) throws Exception{
        
        //basic information, get the zero instance.
        snpName = all_individuals.get(0).getSnpName();
        chromosome = all_individuals.get(0).getChromosome();
        position = all_individuals.get(0).getPosition();
        
       ArrayList<IndividualSnpData> het_individuals = UtilityMethods.isolateValidHeterozygotesFromIndividualSnpData(all_individuals);
        
        
        numberOfHets = het_individuals.size();
        
        hetSampleNames = new ArrayList<String>();
        
        asRef = new ArrayList<Integer>(); 
        asAlt = new ArrayList<Integer>(); 
        asNo  = new ArrayList<Integer>(); 
        
        cellProp = new ArrayList<Double>();
        
                
        int total_overlap = 0;
        int ref_total = 0;
        
        for (IndividualSnpData temp_het : het_individuals) {
            //Do nothing if there is no data in het_individuals

            hetSampleNames.add(temp_het.getSampleName());
            
            asRef.add(temp_het.getRefNum());
            asAlt.add(temp_het.getAltNum());
            asNo.add(temp_het.getNoNum());
            
            cellProp.add(temp_het.getCellTypeProp());
            
            //this is used to check if we will continue with calculations.
            total_overlap += temp_het.getRefNum() + temp_het.getAltNum();
            ref_total += temp_het.getRefNum();
        }
        
        
        if((total_overlap >= GlobalVariables.minReads) && (numberOfHets >= GlobalVariables.minHets) ){
            
            
            
            

            Integer[] asRefArray = asRef.toArray(new Integer[asRef.size()]);
            Integer[] asAltArray = asAlt.toArray(new Integer[asAlt.size()]);
            Double[]  cellPropArray  = cellProp.toArray(new Double[cellProp.size()]);

            // Just using the alternative likelihood of the binomial test here.
            // As out null loglikelihood.

            
            
            try{
                
                CTSnullBinomialLikelihood CTSbinomNull = new CTSnullBinomialLikelihood(asRefArray, asAltArray, cellPropArray);             
               
            
                
                double nullProp = ( (double)ref_total / (double)total_overlap);
                
               
                
                double[] nullInput = { nullProp };
                
                nullLogLik = CTSbinomNull.value(nullInput);

                if(GlobalVariables.verbosity >= 100){
                    System.out.println("LogLik of NULL");
                    System.out.println("\tResidual ratio:                   " + Double.toString(nullInput[0]));
                    System.out.println("\tnullLogLik:                       " + Double.toString(nullLogLik));

                }



                // Now do MLE of the alternative test 
                // more computationally intensive 
                // than the previous one.

                CTSaltBinomialLikelihood CTSbinom = new CTSaltBinomialLikelihood(asRefArray, asAltArray, cellPropArray);
                
                NelderMeadSimplex simplex;
                simplex = new NelderMeadSimplex(2, 1.0, 1.0, 2.0, 0.25, 0.25);
                SimplexOptimizer optimizer = new SimplexOptimizer(GlobalVariables.simplexThreshold, GlobalVariables.simplexThreshold); 
                PointValuePair solutionAlt = optimizer.optimize(
                                                new ObjectiveFunction(CTSbinom),
                                                new MaxEval(GlobalVariables.maximumIterations),
                                                simplex,
                                                GoalType.MINIMIZE,
                                                new InitialGuess(new double[] {0, nullInput[0]}), 
                                                new SearchInterval(0.0, 1.0)
                                                );

                double[] valueAlt = solutionAlt.getPoint();
                
                
                
                
                /*
                    In simulations it was found that in about 30% of the cases
                    A solution was not found, so we do it again if there 
                */
                
                
                if((optimizer.getIterations() <= 5) ){
                    
                    if(GlobalVariables.verbosity >= 100){
                        
                        System.out.println("\nfirst starting point was already in a minimum.");
                        System.out.println("Trying with other starting values.");
                    
                    }
                    
                    ArrayList<double[]> StartingValueList = new ArrayList<double[]>();
                    //These are the starting values to try
                    
                StartingValueList.add(new double[] {0.0  , nullProp});
                StartingValueList.add(new double[] {CTSlinearRegressionResults.slope , CTSlinearRegressionResults.intercept});

//                StartingValueList.add(new InitialGuess(new double[] {0.5  , 0.1       }));
//                StartingValueList.add(new InitialGuess(new double[] {-0.01, 0.7       }));
//                StartingValueList.add(new InitialGuess(new double[] {0.25 , 0.5       }));
                    
                    
                    for(double[] startingValue :StartingValueList ){
                        SimplexOptimizer newOptimizer = new SimplexOptimizer(GlobalVariables.simplexThreshold, GlobalVariables.simplexThreshold); 
                        PointValuePair newSolutionAlt = newOptimizer.optimize(
                                                        new ObjectiveFunction(CTSbinom),
                                                        new MaxEval(GlobalVariables.maximumIterations),
                                                        simplex,
                                                        GoalType.MINIMIZE,
                                                        new InitialGuess(startingValue), 
                                                        new SearchInterval(0.0, 1.0)
                                                        );

                        double[] newValueAlt = newSolutionAlt.getPoint();


                        if(CTSbinom.value(newValueAlt) < CTSbinom.value(valueAlt)) {

                            if(GlobalVariables.verbosity >= 100 && newOptimizer.getIterations() >= 5 ){
                                System.out.println("New starting values are a better fit to the data");
                                System.out.println("keeping results of the new starting values.\n");
                            }
                            //Asign the new values to the actual used ones.
                            valueAlt = newValueAlt;
                            optimizer = newOptimizer; 
                            
                            break;
                        }
                    }
                }
                    
                
                
                altLogLik  = CTSbinom.value(valueAlt);
                iterations = optimizer.getIterations();

                MLEratioCellType = valueAlt[1];
                MLEratioResidual = valueAlt[0];




                //chi squared statistic is determined based on both null and alt loglikelihoods.
                chiSq = LikelihoodFunctions.ChiSqFromLogLik(nullLogLik, altLogLik);

                //determine P value based on distribution
                pVal = LikelihoodFunctions.determinePvalFrom1DFchiSq(chiSq);
                
                if(GlobalVariables.verbosity >= 10){
                    System.out.println("\n--- Starting cell type specific binomial LRT test estimate ---");

                    System.out.println("LogLik of Alt converged to a threshold of " + Double.toString(GlobalVariables.simplexThreshold));
                    System.out.println("\tCelltype ratio:                   " + Double.toString(valueAlt[0]));
                    System.out.println("\tResidual ratio:                   " + Double.toString(valueAlt[1]));
                    System.out.println("\tIterations to converge:           " + Integer.toString(iterations) + "\n");
                    System.out.println("\tNull log likelihood:              " + Double.toString(nullLogLik));   
                    System.out.println("\tAlt log likelihood:               " + Double.toString(altLogLik) + "\n");
                    System.out.println("\tChisq statistic:                  " + Double.toString(chiSq));
                    System.out.println("\tP value:                          " + Double.toString(pVal));
                    System.out.println("--------------------------------------------------------------");
                    
                }
                
                testConverged = true;
                
            } catch(TooManyEvaluationsException e){
                
                if(GlobalVariables.verbosity >= 1){
                    System.out.println("WARNING: Did not converge to a solution for SNP: " + snpName + "in cell type specific beta binomial");
                    System.out.println("         After " + Integer.toString(GlobalVariables.maximumIterations) +   " iterations.");
                    System.out.println("         Continue-ing with the next.");
                    System.out.println("--------------------------------------------------------------");
            
                }
            
            }
            
            // Add some extra values that will make sure that there could be some kind of non-convergence.
            testPerformed = true;
            
        }
        
    
    }
    
    public static String writeHeader(){
        String header = "chr\tpos\tsnpName\tnumHets\tpVal\tchiSq\tcellTypeRatio\tresidualRatio\tnullLogLik\taltLogLik";
        return header;
    }
    
    
    
    public String writeTestStatistics(boolean all_data) {
        
        String out = "";
        
        //snp info
        out += chromosome + "\t";
        out += position + "\t";
        out += snpName + "\t";
        
        if(testPerformed && testConverged){
            out += numberOfHets + "\t";

            // lets have a look at how the decimals behave here, seems all right.
            // otherwise i think I have to do something with printf or something.
            out += Double.toString(pVal) + "\t";
            out += Double.toString(chiSq) + "\t";
            out += Double.toString(MLEratioCellType) + "\t";
            out += Double.toString(MLEratioResidual) + "\t";
            out += Double.toString(nullLogLik) + "\t";
            out += Double.toString(altLogLik);

            


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

    public boolean isTestPerformed() {
        return testPerformed;
    }

}
