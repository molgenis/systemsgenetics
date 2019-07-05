/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

import java.util.ArrayList;
import java.util.Arrays;
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
    private boolean testConverged = false;
    private final boolean outPutAllData = false;
    
    //Test statistics:
    
    //MLE statistics
    int iterations;
    double alphaParam;
    double betaParam;
    
    
    //likelihood of data:
    double nullLogLik;
    double altLogLik;
    double chiSq;
    double pVal;
    
    //precision
    
    double binomRatioCellType;
    double binomRatioResidual;
    int    nulliterations;
    int    altiterations;
    double NullAlphaParam;
    double NullBetaParam;
    double NullbinomRatio;
    
    public CTSBetaBinomialTest(ArrayList<IndividualSnpData> all_individuals, CTSlinearRegression CTSlinearRegressionResults) throws Exception {
    
        //basic information, get the zero instance.
        snpName = all_individuals.get(0).getSnpName();
        chromosome = all_individuals.get(0).getChromosome();
        position = all_individuals.get(0).getPosition();
        
        //isolate heterozygote individuals.
        ArrayList<IndividualSnpData> het_individuals = UtilityMethods.isolateValidHeterozygotesFromIndividualSnpData(all_individuals);
    
        numberOfHets = het_individuals.size();
        
        hetSampleNames = new ArrayList<String>();
        asRef = new ArrayList<Integer>(); 
        asAlt = new ArrayList<Integer>(); 
        asNo  = new ArrayList<Integer>(); 
        
        dispersion = new ArrayList<Double>(); 
        cellProp   = new ArrayList<Double>();
        int total_overlap = 0;
        
        //Get the basic data without doing any tests.
        for (IndividualSnpData temp_het : het_individuals) {
            //Do nothing if there is no data in het_individuals

            hetSampleNames.add(temp_het.getSampleName());
            
            asRef.add(temp_het.getRefNum());
            asAlt.add(temp_het.getAltNum());
            asNo.add(temp_het.getNoNum());
            
            dispersion.add(temp_het.getDispersion());
            cellProp.add(temp_het.getCellTypeProp());
            
            //this is used to check if we will continue with calculations.
            //BASED on the minHets and MinReads
            total_overlap += temp_het.getRefNum() + temp_het.getAltNum();
        }
        
        //Check if we do a test.
        if((total_overlap >= GlobalVariables.minReads) && (numberOfHets >= GlobalVariables.minHets) ){
            // There is data to perform the binomial test, perform it.       

            if(GlobalVariables.verbosity >= 100){
                System.out.println();
              
                System.out.println("SNP name: " + snpName);
                System.out.println("at: chr" + chromosome + ":" + position);
                System.out.println("------------------------------------------------------");
                System.out.println("Num of hets: " + Integer.toString(numberOfHets));
                System.out.println(het_individuals.get(0).getSnpName());
                System.out.println(total_overlap);
                System.out.println("asRef:       " +  asRef.toString());
                System.out.println("asAlt:       " +  asAlt.toString());
                System.out.println("dispersion:  " +  dispersion.toString());
                System.out.println("cellProp:    " +  cellProp.toString());
                System.out.println("Starting non-CTS Beta binomial Null estimation");
            
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
            
            
            
            
            try{
                BetaBinomAltLikelihood betaBinomNull;
                betaBinomNull = new BetaBinomAltLikelihood(asRefArray, 
                                                            asAltArray, 
                                                            dispArray
                                                            );
                NelderMeadSimplex simplex;
                simplex = new NelderMeadSimplex(2);
                SimplexOptimizer optimizer = new SimplexOptimizer(GlobalVariables.simplexThreshold, GlobalVariables.simplexThreshold); //numbers are to which precision you want it to be done.
                PointValuePair solutionNull = optimizer.optimize(
                                                new ObjectiveFunction(betaBinomNull),
                                                new MaxEval(GlobalVariables.maximumIterations),
                                                simplex,
                                                GoalType.MINIMIZE,
                                                new InitialGuess(new double[] {0.5, 0.5}), 
                                                new SearchInterval(0, 1.0)
                                                );

                double[] valueNull = solutionNull.getPoint();

                nullLogLik  = betaBinomNull.value(valueNull);
                nulliterations = optimizer.getIterations();
                NullAlphaParam = valueNull[0];
                NullBetaParam = valueNull[1];
                NullbinomRatio = valueNull[0] / (valueNull[0] + valueNull[1]);

                if(GlobalVariables.verbosity >= 100){
                    System.out.println("LogLik of Null converged to a threshold of " + Double.toString(GlobalVariables.simplexThreshold));
                    System.out.println("\tNull Alpha parameter:      " + Double.toString(valueNull[0]));
                    System.out.println("\tNull Beta parameter:       " + Double.toString(valueNull[1]));
                    System.out.println("\tIterations to converge:    " + Integer.toString(nulliterations) + "\n");
                }

                //CHECK WHAT THE version does in terms of loglik.
                CTSbetaBinomialAltLikelihoodVersion2 CTSbetaBinomAlt;
                CTSbetaBinomAlt = new CTSbetaBinomialAltLikelihoodVersion2(asRefArray, 
                                                        asAltArray,
                                                        dispArray,
                                                        cellPropArray
                                                        );

                // Going to make some initialGuesses, because sometimes 
                // this falls into a local minimum with the MLE, when only the 
                // null parameters are used as a starting point.


                double nonCTSprop = valueNull[0]  / (valueNull[0] + valueNull[1]);


                ArrayList<InitialGuess> GuessList = new ArrayList<InitialGuess>();

                GuessList.add(new InitialGuess(new double[] {0.0  , nonCTSprop}));
                GuessList.add(new InitialGuess(new double[] {CTSlinearRegressionResults.slope , CTSlinearRegressionResults.intercept}));
                GuessList.add(new InitialGuess(new double[] {0.2  , 0.5}));
                GuessList.add(new InitialGuess(new double[] {-0.2 , 0.5}));

                simplex = new NelderMeadSimplex(2);



                ArrayList<double[]> OptimizerResults = new ArrayList<double[]>();
                Double[] OptimizedLogLik  = new Double[8];

                int i=0;
                int lowestIndices = -1;
                double lowestLogLik = 0;


                for(InitialGuess IGuess : GuessList){

                    PointValuePair solutionAlt = optimizer.optimize(
                                            new ObjectiveFunction(CTSbetaBinomAlt),
                                            new MaxEval(500),
                                            simplex,
                                            GoalType.MINIMIZE,
                                            IGuess,
                                            new SearchInterval(0, 1)
                                            );


                    double[] valueAlt = solutionAlt.getPoint();

                    OptimizerResults.add(valueAlt);
                    OptimizedLogLik[i] = CTSbetaBinomAlt.value(valueAlt);


                    //determine lowest loglik.
                    if(i ==0){
                        lowestIndices = 0;
                        lowestLogLik = OptimizedLogLik[i];
                    } else if(OptimizedLogLik[i] < lowestLogLik){
                        lowestIndices = i;
                        lowestLogLik  = OptimizedLogLik[i]; 
                    }

                    if(GlobalVariables.verbosity >= 100){    
                        System.out.println("\nAlt Loglik convergence of starting coordinates" + Arrays.toString(IGuess.getInitialGuess()));
                        System.out.println("\tFinal parameters:     " + Arrays.toString(valueAlt));
                        System.out.println("\tLogLik:               " + Double.toString(OptimizedLogLik[i]));
                    }

                    //only do the last 2 guesses when the celltype ratio is not zero.
                    //This is an accuracy improvement because sometimes a local minimum is found at celtype is 0.0
                    //But it is very unlikely that this happens.
                    if((i >= 2) && (valueAlt[0] != 0.0)){
                        break;
                    }
                    i++;
                    
                }


                //Now select the lowest LogLik.

                double[] bestParams = OptimizerResults.get(lowestIndices);
                altLogLik = OptimizedLogLik[lowestIndices];

                binomRatioCellType = bestParams[0];
                binomRatioResidual = bestParams[1];


                //chi squared statistic is determined based on both null and alt loglikelihoods.
                chiSq = LikelihoodFunctions.ChiSqFromLogLik(nullLogLik, altLogLik);

                //determine P value based on distribution

                pVal = LikelihoodFunctions.determinePvalFrom1DFchiSq(chiSq);

                
                if(GlobalVariables.verbosity >= 10){
                    
                    System.out.println("\n--- Starting cell type specific beta binomial LRT test estimate ---");
                    System.out.println("LogLik of converged to a threshold of " + Double.toString(GlobalVariables.simplexThreshold) + "\n");
                    System.out.println("\tCellType Binomial ratio:       " + Double.toString(binomRatioCellType) + "\n");
                    System.out.println("\tResidual Binomial ratio:       " + Double.toString(binomRatioResidual) + "\n");
                    System.out.println("\tNull log likelihood:           " + Double.toString(nullLogLik));   
                    System.out.println("\tAlt log likelihood:            " + Double.toString(altLogLik) + "\n");
                    System.out.println("\tChisq statistic:               " + Double.toString(chiSq));
                    System.out.println("\tP value:                       " + Double.toString(pVal));
                    System.out.println("-----------------------------------------------------------------------");

                }
                
                 testConverged = true;
           
            } catch(TooManyEvaluationsException e){
                
                if(GlobalVariables.verbosity >= 1){
                    System.out.println("WARNING: Did not converge to a solution for SNP " + snpName + "in cell type specific beta binomial");
                    System.out.println("         After " + Integer.toString(GlobalVariables.maximumIterations) +   " iterations.");
                    System.out.println("         Continue-ing with the next SNP");
                }
            
            }
            
            
            
            //Finally test was done, so we say the test was performed.;
            testPerformed = true;
            
        }
        
    
    }
    
    
    
    /**
     * @return the testPerformed
     */
    public boolean isTestPerformed() {
        return testPerformed;
    }

    public static String writeHeader(){
        String header = "chr\tpos\tsnpName\tnumHets\tpVal\tchiSq\tcellTypeRatio\tresidualRatio\tnullLogLik\taltLogLik";
        return header;
    }
    
    String writeTestStatistics(boolean all_data) {
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
            out += Double.toString(binomRatioCellType) + "\t";
            out += Double.toString(binomRatioResidual) + "\t";
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
    
}
