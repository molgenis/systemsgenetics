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
public class as_all_individual_per_snp {
    
    private int[] as1; //numbers following reference
    private int[] as2; //numbers following alternative
    private int[] as3; //numbers following no allele
    private char[] genotype; // genotypes of individuals.
    private String snp_name;

    /**
     * @return the as1
     */
    public int[] getAs1() {
        return as1;
    }

    /**
     * @param as1 the as1 to set
     */
    public void setAs1(int[] as1) {
        this.as1 = as1;
    }

    /**
     * @return the as2
     */
    public int[] getAs2() {
        return as2;
    }

    /**
     * @param as2 the as2 to set
     */
    public void setAs2(int[] as2) {
        this.as2 = as2;
    }

    /**
     * @return the as3
     */
    public int[] getAs3() {
        return as3;
    }

    /**
     * @param as3 the as3 to set
     */
    public void setAs3(int[] as3) {
        this.as3 = as3;
    }

    /**
     * @return the genotype
     */
    public char[] getGenotype() {
        return genotype;
    }

    /**
     * @param genotype the genotype to set
     */
    public void setGenotype(char[] genotype) {
        this.genotype = genotype;
    }

    /**
     * @return the snp_name
     */
    public String getSnp_name() {
        return snp_name;
    }

    /**
     * @param snp_name the snp_name to set
     */
    public void setSnp_name(String snp_name) {
        this.snp_name = snp_name;
    }
    
    
    
    
    
}
