/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl4;

/**
 *
 * @author harmjan
 */
public class EQTL4 implements Comparable<EQTL4> {

    private double pvalue = 1d;
    private int snp = -1;
    private int probe = -1;
    
    public EQTL4(int snp, int probe, double pvalue){
        this.snp = snp;
        this.probe = probe;
        this.pvalue = pvalue;
    }
    
    public double getPvalue() {
        return pvalue;
    }

    public int compareTo(EQTL4 t) {
        if(t.equals(this)){
            return 0;
        } else if(t.getPvalue() < pvalue){
            return 1;
        } else {
            return -1;
        }
    }
    
    public boolean equals(EQTL4 t){
        if(t.getSNP() == snp && t.getProbe() == probe && t.getPvalue() == pvalue){
            return true;
        } else {
            return false;
        }
    }

    public int getSNP() {
        return snp;
    }

    public int getProbe() {
        return probe;
    }

    void setPvalue(double pvalue) {
        this.pvalue = pvalue;
    }

    void setSNP(int snp) {
        this.snp = snp;
    }

    void setProbe(int probe) {
        this.probe = probe;
    }
    
}
