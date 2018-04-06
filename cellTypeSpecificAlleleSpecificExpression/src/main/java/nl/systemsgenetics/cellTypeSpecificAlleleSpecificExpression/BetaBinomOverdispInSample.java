/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.univariate.SearchInterval;
/**
 *
 * @author adriaan
 */
public class BetaBinomOverdispInSample {
    
    private String sampleName;
    private double[] overdispersion = new double[1];
    
    
    @SuppressWarnings({"empty-statement", "empty-statement"})
    public BetaBinomOverdispInSample(String filename) throws FileNotFoundException, IOException{
        
        sampleName = filename;
        
        if(GlobalVariables.verbosity >= 10){
            System.out.println("\n--- Starting beta binomial dispersion estimate ---");
            System.out.println("AS File: " + sampleName);
            System.out.println("--------------------------------------------------");
        }
        
        ArrayList<IndividualSnpData> allSnps = new ArrayList<IndividualSnpData>();
        
        File file = new File(filename);
        FileReader fr = new FileReader(file);
        BufferedReader br = new BufferedReader(fr);
        String line;
        while((line = br.readLine()) != null){
            allSnps.add(new IndividualSnpData(file.getAbsolutePath(), line));
        }
        br.close();
        fr.close();
        
        ArrayList<IndividualSnpData> hets;
        hets = UtilityMethods.isolateValidHeterozygotesFromIndividualSnpData(allSnps);
        
        int numOfHets = hets.size();
       
        int[] asRef = new int[numOfHets]; 
        int[] asAlt = new int[numOfHets];
        
        int totalRef = 0;
        int totalAlt = 0;
        
        for(int i =0; i< numOfHets; i ++){
            asRef[i] = hets.get(i).getRefNum();
            asAlt[i] = hets.get(i).getAltNum();
            
            totalRef +=asRef[i];
            totalAlt +=asAlt[i];
                    
            
        }
        if(GlobalVariables.verbosity >= 10){
            System.out.println("sample loaded");

            System.out.println("\ttotal_ref = " + totalRef);
            System.out.println("\ttotal_alt = " + totalAlt);
        }
        
        
        BetaBinomLikelihoodForOverdispersion betaBinom = new BetaBinomLikelihoodForOverdispersion(asRef, asAlt);
        NelderMeadSimplex simplex;
        simplex = new NelderMeadSimplex(1);
        SimplexOptimizer optimizer = new SimplexOptimizer(GlobalVariables.simplexThreshold, GlobalVariables.simplexThreshold); //numbers are to which precision you want it to be done.
        PointValuePair solution = optimizer.optimize(
                                            new ObjectiveFunction(betaBinom),
                                            new MaxEval(GlobalVariables.maximumIterations),
                                            simplex,
                                            GoalType.MINIMIZE,
                                            new InitialGuess(new double[] {0.5}), 
                                            new SearchInterval(0.0, 1.0)
                                            );
        
        overdispersion = solution.getFirst();
        if(GlobalVariables.verbosity >= 10){
            System.out.println("Log likelihood converged to a threshold of " + Double.toString(GlobalVariables.simplexThreshold));
            System.out.println("\tDispersion sigma: " + Double.toString(overdispersion[0]));
            System.out.println("\tLog likelihood:   " + Double.toString(betaBinom.value(overdispersion)));
        }

        

    }
    
    
     public BetaBinomOverdispInSample(String sampleName, double[] dispersion ){
         this.sampleName = sampleName;
         this.overdispersion = dispersion;
     }
    
    

    /**
     * @return the sampleName
     */
    public String getSampleName() {
        return sampleName;
    }

    /**
     * @return the overdispersion
     */
    public double[] getOverdispersion() {
        return overdispersion;
    }


    

}
