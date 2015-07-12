/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

/**
 *
 * @author adriaan
 */
class CTSindividualSnpData extends IndividualSnpData{

    private double cellTypeProp;
    
    public CTSindividualSnpData(String SAMPLENAME, String snpLine, double cellTypeProp1) {
        super(SAMPLENAME, snpLine);
        
        cellTypeProp = cellTypeProp1;
        
    }

    /**
     * @return the cellTypeProp
     */
    public double getCellTypeProp() {
        return cellTypeProp;
    }
    
    
    
}
