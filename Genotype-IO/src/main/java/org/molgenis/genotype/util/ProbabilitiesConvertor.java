/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.util;

import java.util.Arrays;
import java.util.List;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.probabilities.SampleVariantProbabilities;

/**
 *
 * @author Patrick Deelen
 */
public class ProbabilitiesConvertor {

	/**
	 * Convert alleles to probability values. If more than two alleles will
	 * return all missing. If sample alleles found not specified by alleles then sample will be set to missing
	 *
	 *
	 * @param sampleAlleles
	 * @param alleles
	 * @return
	 */
	public static SampleVariantProbabilities[] convertCalledAllelesToDosage(List<Alleles> sampleAlleles, Alleles alleles) {

		SampleVariantProbabilities[] probs = new SampleVariantProbabilities[sampleAlleles.size()];

		if (alleles.getAlleleCount() == 0) {
			throw new GenotypeDataException("Error converting alleles to probabilities. No alleles detected");
		} else if (alleles.getAlleleCount() > 2) {

			Arrays.fill(probs, SampleVariantProbabilities.MISSING_PROB);
			return probs;

		} else {

			Alleles AA = Alleles.createAlleles(alleles.get(0), alleles.get(0));
			Alleles AB;
			Alleles BA;
			Alleles BB;

			if (alleles.getAlleleCount() == 2) {
				AB = Alleles.createAlleles(alleles.get(0), alleles.get(1));
				BA = Alleles.createAlleles(alleles.get(1), alleles.get(0));
				BB = Alleles.createAlleles(alleles.get(1), alleles.get(1));
			} else {
				AB = Alleles.createAlleles(Allele.ZERO, Allele.ZERO);
				BA = Alleles.createAlleles(Allele.ZERO, Allele.ZERO);
				BB = Alleles.createAlleles(Allele.ZERO, Allele.ZERO);
			}

			for (int i = 0; i < probs.length; ++i) {
				Alleles sampleVariant = sampleAlleles.get(i);

				if (sampleVariant.contains(Allele.ZERO)) {
					probs[i] = SampleVariantProbabilities.MISSING_PROB;
				} else if (sampleVariant == AA) {
					probs[i] = SampleVariantProbabilities.AA_PROB;
				} else if (sampleVariant == AB) { //We test above for missing so is never true for monomorphic SNPs
					probs[i] = SampleVariantProbabilities.AB_PROB;
				} else if (sampleVariant == BA) { //We test above for missing so is never true for monomorphic SNPs
					probs[i] = SampleVariantProbabilities.AB_PROB;
				} else if (sampleVariant == BB) { //We test above for missing so is never true for monomorphic SNPs
					probs[i] = SampleVariantProbabilities.BB_PROB;
				} else {
					probs[i] = SampleVariantProbabilities.MISSING_PROB;
				}
				
			}
	
			return probs;
			
		}
		
	}
}
