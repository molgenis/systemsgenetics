/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cell_type_specific_ase;

/**
 *
 * @author adriaan
 * Simple class adding the double dispersionData to the IndividualSNP data, for usage in the Beta Binomial Test.
 */
public class dispersedIndividualSnpData extends IndividualSnpData{
    private double dispersion;
    
    public dispersedIndividualSnpData(String SAMPLENAME, String snpLine, double dispersionParameter) {
        super(SAMPLENAME, snpLine);
        dispersion = dispersionParameter;
    }
    
    public double getDispersion(){
        return dispersion;
    }
}
