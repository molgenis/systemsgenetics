/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.special.Beta;
import org.jdom.IllegalDataException;

/**
 *
 * @author adriaan
 */

public class LikelihoodFunctions {
    
    static double BinomLogLik(double d, int asRef, int asAlt) {
        int total = asRef + asAlt;
        //determine likelihood here.
        BinomialDistribution binomDist = new BinomialDistribution(total, d);
        Double logDensity = binomDist.logProbability(asRef);
        return -1.0 * logDensity;
    }
    
    static double BinomLogLik(double d, int[] asRef, int[] asAlt) {
        
        double logLik = 0; 
        
        for(int i = 0; i < asRef.length; i++){
            //determine likelihood here.
           
            logLik += BinomLogLik(d, asRef[i], asAlt[i]);
        }
        return logLik;
    }
    
    
     static double BetaBinomLogLik(double sigma, double alpha, double beta, int asRef, int asAlt){
        //Note, hetp and error are taken from WASP, van de Geijn et al. (2015)
        int AS1      = asRef;
        int AS2      = asAlt;
        double hetp  = GlobalVariables.hetProb;
        double error = GlobalVariables.seqError;

        double a;
        double b;

        a = (alpha) /(alpha + beta) *  (1.0 / Math.pow(sigma, 2) - 1.0);
        b = (beta)  /(alpha + beta) *  (1.0 / Math.pow(sigma, 2) - 1.0);
        
        
        double part1 = 0.0;
        part1 += Beta.logBeta(AS1 + a, AS2 + b);
        part1 -= Beta.logBeta(a, b);
        
        
        if(GlobalVariables.hetProb == 1.0){
            return -1.0 * part1;
        }
        
        //Below is unused in our case. hetProb is always 1.0
        double e1 = Math.log(error) * AS1 + Math.log(1.0 - error) * AS2;
        double e2 = Math.log(error) * AS2 + Math.log(1.0 - error) * AS1;
        
        if(GlobalVariables.hetProb == 0.0){
            return -1.0 * addlogs(e1, e2);
        }
        
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
        if(chiSq < 0){
            //sometimes the value is not really 0 but something in the order 10e-15,
            //probably due to floating point errors.
            //throws an exception if the negative value is far from 0
            if(chiSq < -0.00001){
                //This was done because Freerk encountered this error.
                System.err.println(Double.toString(chiSq));
                throw new IllegalDataException("ChiSq value is lower than 1.");
            }
            chiSq = 0.0;
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
