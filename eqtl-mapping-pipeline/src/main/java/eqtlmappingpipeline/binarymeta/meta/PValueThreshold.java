/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.binarymeta.meta;

/**
 *
 * @author harmjan
 */
public class PValueThreshold {

    private double pvalue = Double.MAX_VALUE;

    /**
     * @return the pvalue
     */
    public double getPvalue() {
	return pvalue;
    }

    /**
     * @param pvalue the pvalue to set
     */
    public void setPvalue(double pvalue) {
	this.pvalue = pvalue;
    }
}
