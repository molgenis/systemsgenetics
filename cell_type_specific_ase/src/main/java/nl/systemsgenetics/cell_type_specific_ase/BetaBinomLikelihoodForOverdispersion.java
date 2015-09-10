/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cell_type_specific_ase;

import org.apache.commons.math3.special.Beta;
import org.ejml.simple.SimpleMatrix;

/**
 * 
 * @author adriaan
 * 
 */

public class BetaBinomLikelihoodForOverdispersion  implements Function  {
    
    private int[] asRef;
    private int[] asAlt;
    
    
    
    public BetaBinomLikelihoodForOverdispersion(int[] asRef1, int[] asAlt1){
        
        asRef = asRef1;
        asAlt = asAlt1;
    
    }
    
    
    @Override
    public int getDim(){return 1;}
    
    
    @Override
    public double value(double[] t) {
        
        if(t.length != 1 ){
            throw new RuntimeException("BetaBinom function requires 1 input!");
        }
        
        if(asRef.length == 0){
            return 0.0;
        }
        
        if(asRef.length != 0){
        
        
        }
        
        
        double logLik = 0.0;
        
        
        
        for(int i=0; i < asRef.length; i++ ){
            double sigma = t[0];
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
    
    @Override
    public double[] value(SimpleMatrix xx){
        
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
    
    
    public double validate_with_output(double[] t) {
        
        if(t.length != 1 ){
            throw new RuntimeException("BetaBinom function requires 1 input!");
        }
        
        if(asRef.length == 0){
            return 0.0;
        }
        
        if(asRef.length != 0){
        
        
        }
        
        
        double logLik = 0.0;
        
        
        
        for(int i=0; i < asRef.length; i++ ){
            double sigma = t[0];
            double AS1   = asRef[i];
            
            double AS2   = asAlt[i];
            double hetp  = 0.980198;
            double error = 0.005;
                    
            double a;
            double b;
            System.out.println("------\nNew iteration");
            a = Math.exp( Math.log(0.5) + Math.log((1.0 / Math.pow(sigma, 2)) - 1.0));
            System.out.println("a: " + Double.toString(a));
            
            b = Math.exp( Math.log(0.5) + Math.log((1.0 / Math.pow(sigma, 2)) - 1.0));
            System.out.println("b: " + Double.toString(b));

            double part1 = 0.0;
            part1 += Beta.logBeta(AS1 + a, AS2 + b);
            System.out.println("Beta.logBeta(" + "AS1 + a, AS2 + b)" + Double.toString(Beta.logBeta(AS1 + a, AS2 + b)));
            part1 -= Beta.logBeta(a, b);
            System.out.println("Beta.logBeta( a,b)" + Double.toString(Beta.logBeta(a, b)));
            System.out.println("part1: " + Double.toString(part1));

            double e1 = Math.log(error) * AS1 + Math.log(1.0 - error) * AS2;
            System.out.println("e1:" + Double.toString(e1));

            double e2 = Math.log(error) * AS2 + Math.log(1.0 - error) * AS1;
            System.out.println("e2:" + Double.toString(e2));

            logLik +=  addlogs(Math.log(hetp) + part1, Math.log(1 - hetp) + addlogs(e1, e2));
            System.out.println("addlogs(Math.log(hetp) + part1, Math.log(1 - hetp) + addlogs(e1, e2)): " + Double.toString(addlogs(Math.log(hetp) + part1, Math.log(1 - hetp) + addlogs(e1, e2))));
            System.out.println("Ended iteration\n------");
        }
        
        return -1.0 * logLik;
    }
    
    
    
}
