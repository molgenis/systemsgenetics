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

            int AS1   = asRef[i];

            int AS2   = asAlt[i];
            
            logLik +=  LikelihoodFunctions.BetaBinomLogLik(sigma, alpha, beta, AS1, AS2);
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
    
    
    
    private double addlogs(double loga, double logb){
        double i = Math.max(loga, logb) + Math.log(1 + Math.exp(-Math.abs(loga -logb)));
        return i;
    
    }
    
}
