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
class CTSdispersedIndividualSnpData extends IndividualSnpData{

    private double cellTypeProp;
    private double dispersion;
    
    public CTSdispersedIndividualSnpData(String SAMPLENAME, String snpLine, double cellTypeProp1, double dispersion1) {
        super(SAMPLENAME, snpLine);
        
        dispersion = dispersion1;
        cellTypeProp = cellTypeProp1;
        
    }

    /**
     * @return the cellTypeProp
     */
    public double getCellTypeProp() {
        return cellTypeProp;
    }
    
    public double getDispersion(){
        return dispersion;
    }
    
    
}
