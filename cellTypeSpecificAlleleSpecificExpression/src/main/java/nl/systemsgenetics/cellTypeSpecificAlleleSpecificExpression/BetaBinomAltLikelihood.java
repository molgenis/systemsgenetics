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
