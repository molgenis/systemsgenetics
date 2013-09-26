/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.probabilities;

/**
 *
 * @author Patrick Deelen
 */
public interface SampleVariantProbabilities {

	public static final SampleVariantProbabilities MISSING_PROB = new SampleVariantProbabilities3Probs(new float[]{0f, 0f, 0f});
	public static final SampleVariantProbabilities AA_PROB = new SampleVariantProbabilities3Probs(new float[]{1f, 0f, 0f});
	public static final SampleVariantProbabilities AB_PROB = new SampleVariantProbabilities3Probs(new float[]{0f, 1f, 0f});
	public static final SampleVariantProbabilities BB_PROB = new SampleVariantProbabilities3Probs(new float[]{0f, 0f, 1f});

	/**
	 * Vector of 3 probabilities AA AB BB
	 *
	 * @return
	 */
	float[] getProbilities();
}
