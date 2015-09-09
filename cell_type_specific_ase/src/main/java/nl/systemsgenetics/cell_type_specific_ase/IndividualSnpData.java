/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cell_type_specific_ase;

/**
 *
 * @author adriaan
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
    
    
   public IndividualSnpData(String SAMPLENAME, String snpLine ){
       
       sampleName = SAMPLENAME;
       String[] values = snpLine.split("\t");
       
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


}
