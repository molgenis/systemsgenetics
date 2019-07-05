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
public class CTSbetaBinomialAltLikelihoodVersion2 implements Function {

    private Integer[] asRef;
    private Integer[] asAlt;
    private Double[]  dispersion;
    private Double[]  cellProp;
    
    public CTSbetaBinomialAltLikelihoodVersion2(Integer[] asRef1, Integer[] asAlt1, Double[] dispersion1, Double[] cellProp1){
        
        asRef = asRef1;
        asAlt = asAlt1;
        dispersion = dispersion1;
        cellProp = cellProp1;
        
    }    
    
    
    @Override
    public int getDim() {
        return 2;
    }
    
    @Override
    public double value(double[] t){
        double logLik = 0.0;
        
        if(t.length != 2 ){
            throw new RuntimeException("Function requires 2 inputs!");
        }
        
        //Make sure the value is within bounds
        double CellTypeProp = t[0];
        double ResidualProp  = t[1];
        
        for(int i=0; i < asRef.length; i++ ){
            
            
            double fullProp = CellTypeProp  * cellProp[i] + ResidualProp;
            
            if(fullProp > 1 || fullProp < 0){
                return Double.MAX_VALUE;
            }
            
            
            
            double sigma = dispersion[i];
            int    AS1   = asRef[i]; 
            int    AS2   = asAlt[i];
            
            /*
             * here I turn fullprop into alpha and beta parameters
             * assuming alpha is 1, i change beta into the other proportion,
             * while keeping the following assumption: 
             *  prop = alpha / (alpha + beta)
             *  alpha = 1
             *  prop = 1 / (1 + beta)
             *  beta = (1/prop) - 1 
             * 
             * 
             * the python script below shows that this is the case.
            
            for i in range(2000):
                prop = random.uniform(0.0,1.0)
                alpha = 1
                beta = (alpha / prop) - alpha
                reflogps = [math.log(alpha) - math.log(alpha + beta),math.log(beta) - math.log(alpha + beta) ]
                for i in range(2000):  
                    alpha = random.uniform(0, 100)
                    beta  = (alpha / prop) - alpha
                    logps = [math.log(alpha) - math.log(alpha + beta),math.log(beta) - math.log(alpha + beta) ]
                    numpy.testing.assert_array_almost_equal(reflogps, logps, decimal=7)

             * 
             * 
             */
            
            logLik +=  LikelihoodFunctions.BetaBinomLogLik(sigma, 1.0, (1.0 / fullProp) -1 , AS1, AS2);
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
