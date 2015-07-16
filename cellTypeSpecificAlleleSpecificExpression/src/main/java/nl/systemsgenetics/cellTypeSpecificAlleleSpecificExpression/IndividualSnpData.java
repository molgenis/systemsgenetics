/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

import org.jdom.IllegalDataException;

/**
 *
 * @author adriaan
 * This class contains all the data of an individual SNP. 
 * It is creating by passing a line of all the files through it.
 * 
 */
public class IndividualSnpData {
    final String sampleName;
    final String snpName;
    final String chromosome;
    final String position;
    
    
    final char   reference;
    final char   alternative;
    
    final int    refNum;
    final int    altNum;
    final int    noNum;

    final String genotype;
    
    //These values are used to check if the thing has dispersion or cellprop
    
    private boolean hasDispersion = false ;
    private boolean hasCellCount = false;
    
    //These values are optional.
    private double cellTypeProp;
    private double dispersion;
    
    
   public IndividualSnpData(String SAMPLENAME, String snpLine ){
       
       sampleName = SAMPLENAME;
       String[] values = snpLine.split("\t");
       
       if(values.length != 9){
           System.out.println(snpLine);
           throw new IllegalDataException("The SNP line printed above did not contain 9 elements.");
       }
       
       chromosome = values[0];
       position   = values[1];
       snpName    = values[2];
       
       reference    = values[3].charAt(0); 
       alternative  = values[4].charAt(0);
       
       refNum = Integer.parseInt(values[5]);
       altNum = Integer.parseInt(values[6]);
       noNum  = Integer.parseInt(values[7]);
   
       genotype = values[8];
   }
   
    /**
     * @param cellTypeProp the cellTypeProp to set
     * Changes the hasCellCount field as well, for testing.
     */
    public void setCellTypeProp(double cellTypeProp) {
        this.cellTypeProp = cellTypeProp;
        hasCellCount = true;
    }
   
    /**
     * @param dispersion the dispersion to set
     * Changes the hasCellCount field as well.
     */
    public void setDispersion(double dispersion) {
        this.dispersion = dispersion;
        hasDispersion = true;
    }

    //created automatically below.

    /**
     * @return the sampleName
     */
    public String getSampleName() {
        return sampleName;
    }

    /**
     * @return the snpName
     */
    public String getSnpName() {
        return snpName;
    }

    /**
     * @return the chromosome
     */
    public String getChromosome() {
        return chromosome;
    }

    /**
     * @return the position
     */
    public String getPosition() {
        return position;
    }

    /**
     * @return the reference
     */
    public char getReference() {
        return reference;
    }

    /**
     * @return the alternative
     */
    public char getAlternative() {
        return alternative;
    }

    /**
     * @return the refNum
     */
    public int getRefNum() {
        return refNum;
    }

    /**
     * @return the altNum
     */
    public int getAltNum() {
        return altNum;
    }

    /**
     * @return the noNum
     */
    public int getNoNum() {
        return noNum;
    }

    /**
     * @return the genotype
     */
    public String getGenotype() {
        return genotype;
    }

    /**
     * @return the cellTypeProp
     */
    public double getCellTypeProp() throws Exception {
        if(hasCellCount){
            return cellTypeProp;
        }else{
            throw new Exception("CellProp was not added in SNP data, you cannot retrieve it");
        
        }
    }

   

    /**
     * @return the dispersion
     */
    public double getDispersion() throws Exception {
        if(hasDispersion){
            return dispersion;
        } else{
            throw new Exception("Dispersion was not added in this SNP data, you cannot retrieve it");
        }
    }


    /**
     * @return the hasDispersion
     */
    public boolean hasDispersion() {
        return hasDispersion;
    }

    /**
     * @return the hasCellCount
     */
    public boolean hasCellCount() {
        return hasCellCount;
    }


}
