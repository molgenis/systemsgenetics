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
public class GenotypeRecordGt implements GenotypeRecord {

	private final Alleles alleles;

	public GenotypeRecordGt(Alleles alleles) {
		this.alleles = alleles;
	}
	
	@Override
	public Object getGenotypeRecordData(String recordId) {
		if(recordId.equals("GT")){
			return alleles;
		} else {
			return null;
		}
	}
	
	@Override
	public Alleles getSampleAlleles() {
		return alleles;
	}

	@Override
	public float[] getSampleProbs() {
		return null;
	}

	@Override
	public float getSampleDosage() {
		return Float.NaN;
	}

	@Override
	public boolean containsGenotypeRecord(String recordId) {
		return recordId.equals("GT");
	}
	
}
