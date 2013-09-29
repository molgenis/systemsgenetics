/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.probabilities;

import org.apache.commons.lang.ArrayUtils;

/**
 *
 * @author Patrick Deelen
 */
public class SampleVariantProbabilities3Probs implements SampleVariantProbabilities{
	
	/**
	 * MUST BE FINAL!!! Do not make mutable. Create new subtype if you want to do this!!!!
	 */
	private final float[] probs;

	public SampleVariantProbabilities3Probs(float[] probs) {
		this.probs = probs;
	}

	@Override
	public float[] getProbilities() {
		return probs;
	}

	@Override
	public SampleVariantProbabilities getReversedProbilities() {
		
		float[] reverse = new float[probs.length];
		
		for(int i = 0 ; i < probs.length ; ++i){
			reverse[probs.length - i -1] = probs[i];
		}
		
		return new SampleVariantProbabilities3Probs(reverse);
	}
	
	
	
}
