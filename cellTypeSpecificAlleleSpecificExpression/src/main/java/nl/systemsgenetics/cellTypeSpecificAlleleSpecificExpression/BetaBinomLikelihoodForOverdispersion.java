/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

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

        
        return LikelihoodFunctions.BetaBinomLogLik(t[0], 0.5, 0.5, asRef, asAlt);
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
    
}
