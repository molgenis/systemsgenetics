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
 */
class BetaBinomNullLikelihood implements Function {

    private Integer[] asRef;
    private Integer[] asAlt;
    private Double[]  dispersion;
    
    
    public BetaBinomNullLikelihood(Integer[] asRef1, Integer[] asAlt1, Double[] dispersion1){
        
        asRef = asRef1;
        asAlt = asAlt1;
        dispersion = dispersion1;
        
    }
    
    
    
    @Override
    public int getDim() {
        return 1;
    }

    @Override
    public double value(double[] t) {
              
        if(t.length != 1 ){
            throw new RuntimeException("BetaBinomNull function requires 1 input!");
        }
        
        if(asRef.length == 0){
            return 0.0;
        }
        
        
        double logLik = 0.0;
        
        for(int i=0; i < asRef.length; i++ ){
            
            double sigma = dispersion[i];


            int AS1   = asRef[i];


            int AS2   = asAlt[i];
            
            logLik +=  LikelihoodFunctions.BetaBinomLogLik(sigma, t[0], t[0], AS1, AS2);
        }
        
        return logLik;
        
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

    
}
