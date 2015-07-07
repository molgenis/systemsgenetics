/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cell_type_specific_ase;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
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
    
    final String sampleName;
    double[] overdispersion = new double[1];
    
    
    @SuppressWarnings({"empty-statement", "empty-statement"})
    public BetaBinomOverdispInSample(String filename) throws FileNotFoundException, IOException{
        
        sampleName = filename;
        
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
        hets = isolateHeterozygotes(allSnps);
        
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
        
        
        System.out.println("total_ref = " + totalRef);
        System.out.println("total_alt = " + totalAlt);
        
        
        BetaBinomLikelihoodForOverdispersion betaBinom = new BetaBinomLikelihoodForOverdispersion(asRef, asAlt);
        NelderMeadSimplex simplex;
        simplex = new NelderMeadSimplex(1);
        SimplexOptimizer optimizer = new SimplexOptimizer(1e-8, 1e-8); //numbers are to which precision you want it to be done.
        PointValuePair solution = optimizer.optimize(
                                            new ObjectiveFunction(betaBinom),
                                            new MaxEval(500),
                                            simplex,
                                            GoalType.MINIMIZE,
                                            new InitialGuess(new double[] {0.5}), 
                                            new SearchInterval(0.0, 1.0)
                                            );
        
        overdispersion = solution.getFirst();
        
        System.out.println("Overdispersion sigma: " + overdispersion[0]);
        System.out.println("the point of the solution:" + Arrays.toString(solution.getPoint()));
        System.out.println(betaBinom.value(overdispersion));

        
//        I used this for testing compared to python.
//        BetaBinomLikelihood betaBinomTest = new BetaBinomLikelihood(new int[] {4,1}, new int[] {1,4});
//        System.out.println(betaBinomTest.validate_with_output(overdispersion));
//        
        double pythonResult[]; 
        pythonResult= new double[1] ;
        pythonResult[0]= 0.23195;
        System.out.println("python input: " +  betaBinom.value(pythonResult));
        
        //System.out.println("Overdispersion alpha: " + overdispersion[1]);
        //System.out.println("Overdispersion beta: " + overdispersion[2]);

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

    

}
