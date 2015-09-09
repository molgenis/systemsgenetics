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
