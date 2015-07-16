/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

import static java.lang.Math.log;
import java.util.ArrayList;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.special.Beta;

/**
 *
 * @author adriaan
 */
public class likelihoodFunctions {
    
    
    static double BinomLogLik(double d, int[] asRef, int[] asAlt) {
        
        double logLik = 0; 
        
        for(int i = 0; i < asRef.length; i++){
            
            int total = asRef[i]+ asAlt[i];
            //determine likelihood here.
            BinomialDistribution binomDist = new BinomialDistribution(total, d);
            logLik += log(binomDist.probability(asRef[i]));
        
        }
        
        return -1.0 * logLik;
        
    }
    
    
    static double BetaBinomLogLik(double sigma, double alpha, double beta, int[] asRef, int[] asAlt){
        
        double logLik = 0; 
    
        for(int i=0; i < asRef.length; i++ ){
            double AS1   = asRef[i];
            double AS2   = asAlt[i];
            double hetp  = 0.980198;
            double error = 0.005;
                    
            double a;
            double b;

            a = Math.exp( Math.log(0.5) + Math.log((1.0 / Math.pow(sigma, 2)) - 1.0));
            b = Math.exp( Math.log(0.5) + Math.log((1.0 / Math.pow(sigma, 2)) - 1.0));

            double part1 = 0.0;
            part1 += Beta.logBeta(AS1 + a, AS2 + b);
            part1 -= Beta.logBeta(a, b);
            
            double e1 = Math.log(error) * AS1 + Math.log(1.0 - error) * AS2;
            double e2 = Math.log(error) * AS2 + Math.log(1.0 - error) * AS1;
            
            logLik +=  addlogs(Math.log(hetp) + part1, Math.log(1 - hetp) + addlogs(e1, e2));
        }
        
    
        return -1.0 * logLik;
    }
    
    static double addlogs(double loga, double logb){
        double i = Math.max(loga, logb) + Math.log(1 + Math.exp(-Math.abs(loga -logb)));
        return i;
    
    }
    
    static double ChiSqFromLogLik(double nullLogLik, double altLogLik){
        
        double chiSq = 2.0 * (nullLogLik - altLogLik);
        
        //chi sq distribution below cannot handle infinity values.
        //However we just give it the 
        if(chiSq == Double.POSITIVE_INFINITY){
            chiSq = Double.MAX_VALUE;
        }
        
        return chiSq;
    
    }
    
    static double determinePvalFrom1DFchiSq(double chiSq){
    
        //determine P value based on distribution
        ChiSquaredDistribution distribution = new ChiSquaredDistribution(1);
        double pVal = 1 - distribution.cumulativeProbability(chiSq);
        
        return pVal;
    
    }
    
    
    
}
