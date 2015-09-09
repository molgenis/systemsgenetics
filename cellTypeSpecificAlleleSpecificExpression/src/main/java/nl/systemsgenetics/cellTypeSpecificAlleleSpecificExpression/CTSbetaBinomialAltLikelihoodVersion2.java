/*
 * Copyright (C) 2015 Adriaan van der Graaf
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
        
        //Make sure the value is always positive
        double CellTypeProp = t[0];
        double ResidualProp  = t[1];
        
        if(CellTypeProp > 1 || CellTypeProp < 0 ){
            return Double.MAX_VALUE;
        }
        if(ResidualProp > 1 || ResidualProp < 0 ){
            return Double.MAX_VALUE;
        }

        
        // Just some check to make sure everythign is good, but the math is not
        // right at the moment, so I need to check if I can change it.
      
        
        
        
        for(int i=0; i < asRef.length; i++ ){
            
            
            //JUST SOMETHING I'm Trying!! VERY MUCH SOMETHING NOT FINAL!!
            //MATH IS PROPABLY NOT CORRECT.
            //NEED TO CHECK.
            
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
