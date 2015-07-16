/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

import java.util.ArrayList;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
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
    
    //settings for MLE determination
    private final double precision = 1.0E-6;  
    
    
    
    //some settings
    boolean outPutAllData = false;
    boolean testPerformed = false;
    
    public CTSbinomialTest(ArrayList<IndividualSnpData> all_individuals) throws Exception{
        boolean debug=true;
        //basic information, get the zero instance.
        snpName = all_individuals.get(0).getSnpName();
        chromosome = all_individuals.get(0).getChromosome();
        position = all_individuals.get(0).getPosition();
        
       ArrayList<IndividualSnpData> het_individuals = UtilityMethods.isolateHeterozygotesFromIndividualSnpData(all_individuals);
        
        
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
            
            System.out.println();
            System.out.println("--- Starting cell type specific binomial LRT test estimate ---");
            System.out.println("\tSNP name: " + snpName);
            System.out.println("\tat: chr" + chromosome + ":" + position);
            System.out.println("--------------------------------------------------------------");
            
            if(debug){
                
                System.out.println("debug:");
                System.out.println("Num of hets: " + Integer.toString(numberOfHets));
                System.out.println(het_individuals.get(0).getSnpName());
                System.out.println(total_overlap);
                
                System.out.println("asRef:       " +  asRef.toString());
                System.out.println("asAlt:       " +  asAlt.toString());
                System.out.println("cellProp:    " +  cellProp.toString());
            
            }
            
            
            System.out.println("Starting Null estimation");

            Integer[] asRefArray = asRef.toArray(new Integer[asRef.size()]);
            Integer[] asAltArray = asAlt.toArray(new Integer[asAlt.size()]);
            Double[]  cellPropArray  = cellProp.toArray(new Double[cellProp.size()]);

            // Just using the alternative likelihood of the binomial test here.
            // As out null loglikelihood.

            CTSnullBinomialLikelihood CTSbinomNull = new CTSnullBinomialLikelihood(asRefArray, asAltArray, cellPropArray);             
            
            try{
                NelderMeadSimplex simplex;
                simplex = new NelderMeadSimplex(1, 1.0, 1.0, 2.0, 0.25, 0.25);
                SimplexOptimizer optimizer = new SimplexOptimizer(precision, precision); //numbers are to which precision you want it to be done.


                PointValuePair solutionNull = optimizer.optimize(
                                                new ObjectiveFunction(CTSbinomNull),
                                                new MaxEval(20000),
                                                simplex,
                                                GoalType.MINIMIZE,
                                                new InitialGuess(new double[] {1.0 * (ref_total / total_overlap)}),
                                                new SearchInterval(0.0, 1.0)
                                                );

                double[] valueNull = solutionNull.getPoint();
                nullLogLik = CTSbinomNull.value(valueNull);

                if (optimizer.getIterations() < 5) {
                    //I'm assuming that when the number of iterations is too small
                    //The simplex has too big steps to find a solution, so changing the
                    //simplex to make smaller steps. 
                    //Do the optimization again.

                    System.out.println("The number of iterations used to converge was lower than 5.\n"
                                     + "Trying again, with different Simplex.");


                    System.out.println("Old LogLik of nulleNULL converged to a threshold of " + Double.toString(precision));
                    System.out.println("\tResidual ratio:                   " + Double.toString(valueNull[0]));
                    System.out.println("\tIterations to converge:           " + Integer.toString(optimizer.getIterations()) + "\n");


                    simplex = new NelderMeadSimplex(1, 1.0, 1.0, 2.0, 0.1, 0.1);
                    solutionNull = optimizer.optimize(
                                    new ObjectiveFunction(CTSbinomNull),
                                    new MaxEval(20000),
                                    simplex,
                                    GoalType.MINIMIZE,
                                    new InitialGuess(new double[] {valueNull[0]}),
                                    new SearchInterval(0.0, 1.0)
                                    );

                    valueNull = solutionNull.getPoint();
                    nullLogLik = CTSbinomNull.value(valueNull);  
                }


                System.out.println("LogLik of NULL converged to a threshold of " + Double.toString(precision));
                System.out.println("\tResidual ratio:                   " + Double.toString(valueNull[0]));
                System.out.println("\tIterations to converge:           " + Integer.toString(optimizer.getIterations()) + "\n");




                // Now do MLE of the alternative test 
                // Probably more computationally intensive 
                // than the previous test.

                CTSaltBinomialLikelihood CTSbinom = new CTSaltBinomialLikelihood(asRefArray, asAltArray, cellPropArray);

                simplex = new NelderMeadSimplex(2, 1.0, 1.0, 2.0, 0.25, 0.25);
                PointValuePair solutionAlt = optimizer.optimize(
                                                new ObjectiveFunction(CTSbinom),
                                                new MaxEval(500),
                                                simplex,
                                                GoalType.MINIMIZE,
                                                new InitialGuess(new double[] {0, valueNull[0]}), 
                                                new SearchInterval(0.0, 2.0)
                                                );

                double[] valueAlt = solutionAlt.getPoint();

                altLogLik  = CTSbinom.value(valueAlt);
                iterations = optimizer.getIterations();

                MLEratioCellType = valueAlt[0];
                MLEratioResidual = valueAlt[0];




                //chi squared statistic is determined based on both null and alt loglikelihoods.
                chiSq = 2.0 * (nullLogLik - altLogLik);

                //determine P value based on distribution
                ChiSquaredDistribution distribution = new ChiSquaredDistribution(1);
                pVal = 1 - distribution.cumulativeProbability(chiSq);

                System.out.println("LogLik of Alt converged to a threshold of " + Double.toString(precision));
                System.out.println("\tCelltype ratio:                   " + Double.toString(valueAlt[0]));
                System.out.println("\tResidual ratio:                   " + Double.toString(valueAlt[1]));
                System.out.println("\tIterations to converge:           " + Integer.toString(iterations) + "\n");
                System.out.println("\tNull log likelihood:              " + Double.toString(nullLogLik));   
                System.out.println("\tAlt log likelihood:               " + Double.toString(altLogLik) + "\n");
                System.out.println("\tChisq statistic:                  " + Double.toString(chiSq));
                System.out.println("\tP value:                          " + Double.toString(pVal));
                //TODO, I want to format this properly, but not necessary
                System.out.println("\n---- Finished SNP " + snpName);
            
            } catch(TooManyEvaluationsException e){
                System.out.println("WARNING: Did not converge to a solution for " + snpName);
                System.out.println("         After 20,000 iterations.");
                System.out.println("         Continue-ing with the next.");
            }
            
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
            out += Double.toString(MLEratioCellType) + "\t";
            out += Double.toString(MLEratioResidual) + "\t";
            out += Double.toString(nullLogLik) + "\t";
            out += Double.toString(altLogLik) + "\t";

            
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
                    disp_string += Double.toString(cellProp.get(i)) + ";";
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

    public boolean isTestPerformed() {
        return testPerformed;
    }

}
