/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.probabilities;

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
	
}
