/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.variant;

import org.molgenis.genotype.Alleles;

/**
 *
 * @author Patrick Deelen
 */
public class GenotypeRecordProb implements GenotypeRecord {
	
	private final float[] probs;

	public GenotypeRecordProb(float[] probs) {
		this.probs = probs;
	}
	
	@Override
	public Object getGenotypeRecordData(String recordId) {
		if(recordId.equals("GP")){
			return probs;
		} else {
			return null;
		}
	}

	@Override
	public Alleles getSampleAlleles() {
		return null;
	}

	@Override
	public float[] getSampleProbs() {
		return probs;
	}

	@Override
	public float getSampleDosage() {
		return Float.NaN;
	}
	
}
