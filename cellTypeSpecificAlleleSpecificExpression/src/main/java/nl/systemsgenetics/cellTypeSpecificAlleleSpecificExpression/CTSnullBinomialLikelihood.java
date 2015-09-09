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
