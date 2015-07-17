/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

import static java.lang.Math.log;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.special.Beta;

/**
 *
 * @author adriaan
 */
public class likelihoodFunctions {
    
    static double BinomLogLik(double d, int asRef, int asAlt) {
        int total = asRef + asAlt;
        //determine likelihood here.
        BinomialDistribution binomDist = new BinomialDistribution(total, d);
        return -1.0 * log(binomDist.probability(asRef));
    }
    
    static double BinomLogLik(double d, int[] asRef, int[] asAlt) {
        
        double logLik = 0; 
        
        for(int i = 0; i < asRef.length; i++){
            
            int total = asRef[i]+ asAlt[i];
            //determine likelihood here.
           
            logLik += BinomLogLik(d, asRef[i], asAlt[i]);
        }
        return logLik;
    }
    
    
     static double BetaBinomLogLik(double sigma, double alpha, double beta, int asRef, int asAlt){

        int AS1   = asRef;
        int AS2   = asAlt;
        double hetp  = GlobalVariables.hetProb;
        double error = GlobalVariables.seqError;

        double a;
        double b;

        a = Math.exp( (Math.log(alpha) - Math.log(alpha + beta))  + Math.log((1.0 / Math.pow(sigma, 2)) - 1.0));
        b = Math.exp( (Math.log(beta) - Math.log(alpha + beta)) + Math.log((1.0 / Math.pow(sigma, 2)) - 1.0));

        double part1 = 0.0;
        part1 += Beta.logBeta(AS1 + a, AS2 + b);
        part1 -= Beta.logBeta(a, b);

        double e1 = Math.log(error) * AS1 + Math.log(1.0 - error) * AS2;
        double e2 = Math.log(error) * AS2 + Math.log(1.0 - error) * AS1;

        return -1.0 * addlogs(Math.log(hetp) + part1, Math.log(1 - hetp) + addlogs(e1, e2));

     }
    
    static double BetaBinomLogLik(double sigma, double alpha, double beta, int[] asRef, int[] asAlt){
        
        double logLik = 0; 
    
        for(int i=0; i < asRef.length; i++ ){
            int AS1   = asRef[i];
            int AS2   = asAlt[i];
           
            logLik +=  BetaBinomLogLik(sigma, alpha, beta, AS1, AS2);
        }
        
    
        return logLik;
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
