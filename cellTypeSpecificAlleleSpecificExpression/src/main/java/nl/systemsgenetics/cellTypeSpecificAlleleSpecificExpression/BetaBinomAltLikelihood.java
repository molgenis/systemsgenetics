/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

import java.util.ArrayList;
import org.apache.commons.math3.special.Beta;
import org.ejml.simple.SimpleMatrix;

/**
 *
 * @author adriaan
 */
class BetaBinomAltLikelihood implements Function {

    private Integer[] asRef;
    private Integer[] asAlt;
    private Double[] dispersion;
    
    
    public BetaBinomAltLikelihood(Integer[] asRef1, Integer[] asAlt1, Double[] dispersion1){
        
        asRef = asRef1;
        asAlt = asAlt1;
        dispersion = dispersion1;
        
    }
    
    
    
    @Override
    public int getDim() {
        return 2;
    }

    @Override
    public double value(double[] t) {
              
        if(t.length != 2 ){
            throw new RuntimeException("BetaBinomAlt function requires 2 inputs!");
        }
        
        if(asRef.length == 0){
            return 0.0;
        }
        
        
        double logLik = 0.0;
        
        
        double alpha = t[0];
        double beta  = t[1];  
        
        
        for(int i=0; i < asRef.length; i++ ){
            
            double sigma = dispersion[i];
            double AS1   = asRef[i];
            double AS2   = asAlt[i];
            
            //Copied from WASP data.
            double hetp  = 0.980198;
            double error = 0.005;
                    
            double a;
            double b;

            a = Math.exp( (Math.log(alpha) - Math.log(alpha + beta)) + Math.log((1.0 / Math.pow(sigma, 2)) - 1.0));
            b = Math.exp( (Math.log(beta ) - Math.log(alpha + beta)) + Math.log((1.0 / Math.pow(sigma, 2)) - 1.0));

            double part1 = 0.0;
            part1 += Beta.logBeta(AS1 + a, AS2 + b);
            part1 -= Beta.logBeta(a, b);
            
            // If we do not integrate heterozygote calling error or sequencing sequencing error,
            // part1 can be returned in the loop.
            // But now I want to make sure that I follow the WASP approach on the CEU data.
            
            double e1 = Math.log(error) * AS1 + Math.log(1.0 - error) * AS2;
            double e2 = Math.log(error) * AS2 + Math.log(1.0 - error) * AS1;
            
            logLik +=  addlogs(Math.log(hetp) + part1, Math.log(1 - hetp) + addlogs(e1, e2));
        }
        
        return -1.0 * logLik;
        
        
    }

    @Override
    public double[] value(SimpleMatrix xx) {
        int n = xx.numRows();
        double [] retval = new double[n];
        for (int i = 0; i < n; i++){
                retval[i] = value(xx.extractVector(true, i).getMatrix().getData());
        }
		
        return retval;    
    }
    
    
    
    private double addlogs(double loga, double logb){
        double i = Math.max(loga, logb) + Math.log(1 + Math.exp(-Math.abs(loga -logb)));
        return i;
    
    }
    
}
