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
public interface GenotypeRecord {
	
	Object getGenotypeRecordData(String recordId);
	Alleles getSampleAlleles();
	float[] getSampleProbs();
	float getSampleDosage();
	
}
