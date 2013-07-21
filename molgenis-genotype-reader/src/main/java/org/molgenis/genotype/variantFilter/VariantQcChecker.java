/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.variantFilter;

import org.molgenis.genotype.variant.GeneticVariant;

/**
 *
 * @author Patrick Deelen
 */
public class VariantQcChecker {
	
	private double maf;
	private double callRate;
	private double hwe;

	public VariantQcChecker(float maf, float callRate, double hwe) {
		this.maf = maf;
		this.callRate = callRate;
		this.hwe = hwe;
	}
	
	public void setMafCutoff(float maf){
		this.maf = maf;
	}

	public void setHweCutoff(double hwe) {
		this.hwe = hwe;
	}

	public void setCallRateCutoff(double callRate) {
		this.callRate = callRate;
	}
	
	public boolean doesVairantPass(GeneticVariant variant){
		
		if(variant.getMinorAlleleFrequency() < maf){
			return false;
		}
	
		//TODO implement
//		if(variant.getCallRate() < callRate){
//			return false;
//		}
//		
//		if(variant.getHwePvalue() < hwe){
//			return false;
//		}
		
		return true;
		
	}
}
