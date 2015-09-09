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
public class CTSnullBinomialLikelihood implements Function {

    private Integer[] asRef;
    private Integer[] asAlt;
    private Double[]  cellProp;
    
    public CTSnullBinomialLikelihood(Integer[] asRef1, Integer[] asAlt1, Double[] cellProp1){
        asRef = asRef1;
        asAlt = asAlt1;
        cellProp = cellProp1;
    }    
    
    
    @Override
    public int getDim() {
        return 1;
    }
    
    @Override
    public double value(double[] t){
        
        double logLik = 0.0;
        
        if(t.length != 1 ){
            throw new RuntimeException("Function requires 1 inputs!");
        }
        
        double pCellType =  0.0;
        double pResidual = t[0];
        
        double prop;
        
        for(int i=0; i < asRef.length; i++ ){
            
            prop = pCellType * cellProp[i] + pResidual;
            
            // IS THIS CORRECT?
            
            if(prop >= 1){
                prop=1.0;
            }
            if(prop < 0){
                prop=0.0; 
            }
            
            //determine likelihood here.
            logLik += LikelihoodFunctions.BinomLogLik(prop, asRef[i], asAlt[i]);
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
