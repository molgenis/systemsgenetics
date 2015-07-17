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
    double alphaParamCellType;
    double binomRatioCellType;
    double betaParamCellType;
    double alphaParamResidual;
    double betaParamResidual;
    double binomRatioResidual;
    int    nulliterations;
    int    altiterations;
    double NullAlphaParam;
    double NullBetaParam;
    double NullbinomRatio;
    
    public CTSBetaBinomialTest(ArrayList<IndividualSnpData> all_individuals) throws Exception {
    
        boolean debug=true;
        //basic information, get the zero instance.
        snpName = all_individuals.get(0).getSnpName();
        chromosome = all_individuals.get(0).getChromosome();
        position = all_individuals.get(0).getPosition();
        
        //isolate heterozygote individuals.
        ArrayList<IndividualSnpData> het_individuals = UtilityMethods.isolateHeterozygotesFromIndividualSnpData(all_individuals);
    
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
                System.out.println("cellProp:    " +  cellProp.toString());

            
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
            SimplexOptimizer optimizer = new SimplexOptimizer(GlobalVariables.simplexThreshold, GlobalVariables.simplexThreshold); //numbers are to which precision you want it to be done.
            PointValuePair solutionNull = optimizer.optimize(
                                            new ObjectiveFunction(betaBinomNull),
                                            new MaxEval(500),
                                            simplex,
                                            GoalType.MINIMIZE,
                                            new InitialGuess(new double[] {0.5, 0.5}), 
                                            new SearchInterval(-1000.0, 1000.0)
                                            );
            
            double[] valueNull = solutionNull.getPoint();
            
            nullLogLik  = betaBinomNull.value(valueNull);
            nulliterations = optimizer.getIterations();
            NullAlphaParam = valueNull[0];
            NullBetaParam = valueNull[1];
            NullbinomRatio = valueNull[0] / (valueNull[0] + valueNull[1]);
            
            System.out.println("LogLik of Null converged to a threshold of " + Double.toString(GlobalVariables.simplexThreshold));
            System.out.println("\tNull Alpha parameter:      " + Double.toString(valueNull[0]));
            System.out.println("\tNull Beta parameter:       " + Double.toString(valueNull[1]));
            System.out.println("\tIterations to converge:    " + Integer.toString(nulliterations) + "\n");
            
            
            

            
 
            //CHECK WHAT THE version2 DOES in terms of loglik.
            CTSbetaBinomialAltLikelihoodVersion2 CTSbetaBinomAlt;
            CTSbetaBinomAlt = new CTSbetaBinomialAltLikelihoodVersion2(asRefArray, 
                                                    asAltArray,
                                                    dispArray,
                                                    cellPropArray
                                                    );

            simplex = new NelderMeadSimplex(2);
            PointValuePair solutionAlt = optimizer.optimize(
                                            new ObjectiveFunction(CTSbetaBinomAlt),
                                            new MaxEval(500),
                                            simplex,
                                            GoalType.MINIMIZE,
                                            new InitialGuess(new double[] {0.000, valueNull[0]  / (valueNull[0] + valueNull[1]) }), //Start with the loglik of the null.
                                            new SearchInterval(0, 1)
                                            );

            double[] valueAlt = solutionAlt.getPoint();


            altLogLik  = CTSbetaBinomAlt.value(valueAlt);


            altiterations = optimizer.getIterations();


            binomRatioCellType = valueAlt[0];
            binomRatioResidual = valueAlt[1];


            //chi squared statistic is determined based on both null and alt loglikelihoods.
            chiSq = 2.0 * (nullLogLik - altLogLik);

            //determine P value based on distribution
            ChiSquaredDistribution distribution = new ChiSquaredDistribution(1);
            pVal = 1 - distribution.cumulativeProbability(chiSq);

            //remove this line, doesn't always hold.
            System.out.println("LogLik of Alt (version2) converged to a threshold of " + Double.toString(GlobalVariables.simplexThreshold) + "\n");
            System.out.println("\tCellType Binomial ratio:       " + Double.toString(binomRatioCellType) + "\n");
            System.out.println("\tResidual Binomial ratio:       " + Double.toString(binomRatioResidual) + "\n");
            System.out.println("\tIterations to converge:        " + Integer.toString(altiterations) + "\n");
            System.out.println("\tNull log likelihood:           " + Double.toString(nullLogLik));   
            System.out.println("\tAlt log likelihood:            " + Double.toString(altLogLik) + "\n");
            System.out.println("\tChisq statistic:               " + Double.toString(chiSq));
            System.out.println("\tP value:                       " + Double.toString(pVal));
            //TODO, I want to format this properly, but not necessary
            System.out.println("\n---- Finished SNP " + snpName);

 
            
            
            //Finally test was done, so we say the test was performed.
            testPerformed = true;
            
        }
        
    
    }
    
    
    
    /**
     * @return the testPerformed
     */
    public boolean isTestPerformed() {
        return testPerformed;
    }

    String writeTestStatistics(boolean all_data) {
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
            out += Double.toString(alphaParamCellType) + "\t";
            out += Double.toString(betaParamCellType) + "\t";
            out += Double.toString(alphaParamResidual) + "\t";
            out += Double.toString(betaParamResidual) + "\t";
            
            
            if(outPutAllData){
                String samples_string="";
                String ref_string="";
                String alt_string="";                
                String no_string="";
                String disp_string="";
                String cellProp_string="";

                for(int i=0; i < hetSampleNames.size(); i++){
                    
                    //samples_string += hetSampleNames.get(i) + ";";
                    ref_string += Integer.toString(asRef.get(i)) + ";";
                    alt_string += Integer.toString(asAlt.get(i)) + ";";
                    no_string += Integer.toString(asNo.get(i)) + ";";
                    disp_string += Double.toString(dispersion.get(i)) + ";";
                    cellProp_string += Double.toString(cellProp.get(i)) + ";";
                }
                
                //remove last delimiter
                //samples_string = samples_string.substring(0, samples_string.length()-1);
                ref_string = ref_string.substring(0, ref_string.length()-1);
                alt_string = alt_string.substring(0, alt_string.length()-1);
                no_string = no_string.substring(0, no_string.length()-1);
                disp_string = disp_string.substring(0, disp_string.length()-1);
                cellProp_string = cellProp_string.substring(0, cellProp_string.length()-1);
                
                //out += "\t" + samples_string + "\t" + ref_string + "\t" + alt_string + "\t" + no_string;
                out += "\t" + ref_string + "\t" + alt_string + "\t" + no_string  + "\t" + disp_string + "\t" + cellProp_string;


            }


        } else {
            //when no testing is done, will only output snp name and position, and NA.
            //Make sure this is still correct.
            
            for(int i=0; i < 8; i++ ){
                out += "NA\t";
            
            }
            if(outPutAllData){
                for(int i=0; i < 6; i++ ){
                    out += "NA\t";

                }
            }
            out += "NA";
        
        }
        
        return out;   
        
    }
    
}
