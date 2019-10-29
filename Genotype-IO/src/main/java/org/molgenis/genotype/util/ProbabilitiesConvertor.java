/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;

import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;

/**
 *
 * @author Patrick Deelen
 */
public class ProbabilitiesConvertor {

	/**
	 * Convert alleles to probability values. If more than two alleles will
	 * return all missing. If sample alleles found not specified by alleles then
	 * sample will be set to missing
	 *
	 *
	 * @param sampleAlleles
	 * @param alleles
	 * @return
	 */
	public static float[][] convertCalledAllelesToProbability(List<Alleles> sampleAlleles, Alleles alleles) {

		float[][] probs = new float[sampleAlleles.size()][3];

		if (alleles.getAlleleCount() == 0) {
			throw new GenotypeDataException("Error converting alleles to probabilities. No alleles detected");
		} else if (alleles.getAlleleCount() > 2) {

			Arrays.fill(probs, new float[]{0,0,0});
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
					probs[i] = new float[]{0,0,0};
				} else if (sampleVariant == AA) {
					probs[i] = new float[]{1,0,0};
				} else if (sampleVariant == AB) { //We test above for missing so is never true for monomorphic SNPs
					probs[i] = new float[]{0,1,0};
				} else if (sampleVariant == BA) { //We test above for missing so is never true for monomorphic SNPs
					probs[i] = new float[]{0,1,0};
				} else if (sampleVariant == BB) { //We test above for missing so is never true for monomorphic SNPs
					probs[i] = new float[]{0,0,1};
				} else {
					probs[i] = new float[]{0,0,0};
				}

			}

			return probs;

		}

	}

	/**
	 * Uses inexact heuristic as defined by plinkseq
	 * (thtp://atgu.mgh.harvard.edu/plinkseq/dosage.shtml) to convert dosage to
	 * probabilities.
	 *
	 * @param sampleDosages assuming count of ref allele.
	 * @return
	 */
	public static float[][] convertDosageToProbabilityHeuristic(float[] sampleDosages) {

		float[][] probs = new float[sampleDosages.length][3];

		for (int i = 0; i < sampleDosages.length; ++i) {

			float sampleDosage = sampleDosages[i];

			if (sampleDosage > 2 || sampleDosage < 0) {
				probs[i] = new float[]{0f, 0f, 0f};
			} else if (sampleDosage < 1) {
				probs[i] = new float[]{0, sampleDosage, 1 - sampleDosage};
			} else {
				//sampleDosage >= 1 && sampleDosage <= 2
				probs[i] = new float[]{sampleDosage - 1, 2 - sampleDosage, 0};
			}

		}

		return probs;

	}

	/**
	 *
	 * @param probs
	 * @param variantAlleles the two alleles for this variant
	 * @param minProbability to call a genotype
	 * @return
	 */
	public static List<Alleles> convertProbabilitiesToAlleles(float[][] probs, Alleles variantAlleles, double minProbability) {

		ArrayList<Alleles> sampleAlleles = new ArrayList<Alleles>(probs.length);

		final int alleleCount = variantAlleles.getAlleleCount();

		if (alleleCount > 2 || alleleCount == 0) {
			throw new GenotypeDataException("Error converting posterior probabilities to called alleles. Found non biallelic SNP");
		}

		Alleles aa = Alleles.createAlleles(variantAlleles.get(0), variantAlleles.get(0));
		Alleles bb;

		if (alleleCount == 2) {
			bb = Alleles.createAlleles(variantAlleles.get(1), variantAlleles.get(1));
		} else {
			bb = null;
		}

		Alleles missing = Alleles.createAlleles(Allele.ZERO, Allele.ZERO);

		for (float[] sampleProbs : probs) {

			int maxProbIndex = -1;
			float maxProb = 0;

			int i = 0;
			for (float prob : sampleProbs) {
				if (prob > 0 && prob >= minProbability && prob > maxProb) {
					maxProbIndex = i;
					maxProb = prob;
				}
				++i;
			}

			if (alleleCount == 1 && maxProbIndex >= 1) {
				throw new GenotypeDataException("Error converting posterior probabilities to called alleles. Illigale probability.");
			}

			switch (maxProbIndex) {
				case -1:
					sampleAlleles.add(missing);
					break;
				case 0:
					sampleAlleles.add(aa);
					break;
				case 1:
					sampleAlleles.add(variantAlleles);
					break;
				case 2:
					sampleAlleles.add(bb);
					break;
				default:
					throw new GenotypeDataException("Error converting posterior probabilities to called alleles. This should not happen, please report this bug.");
			}

		}

		return sampleAlleles;

	}

	public static float[] convertProbabilitiesToDosage(float[][] probs, double minProbability) {

		float[] dosages = new float[probs.length];

		for (int i = 0; i < probs.length; ++i) {

			boolean containsMinProbability = false;

			for (float prob : probs[i]) {
				if (prob >= minProbability) {
					containsMinProbability = true;
					break;
				}
			}

			if (containsMinProbability) {

				dosages[i] = (probs[i][0] * 2) + probs[i][1];
				if (dosages[i] > 2) {
					dosages[i] = 2;
				}

			} else {
				dosages[i] = -1;
			}
		}

		return dosages;


	}

	/**
	 * Method for converting BGEN probabilities to regular posterior probabilities.
	 * unlike the regular posterior probabilities, probabilities in BGEN can be stored
	 * in an array of arbitrary length corresponding to the number of ordered combinations of
	 * alleles for a given ploidy. Any number of alleles can be represented.
	 *
	 * If the BGEN probabilities for a sample represent 3 possible genotypes (biallelic for dipoid samples)
	 * This is returned.
	 *
	 * If the BGEN probabilities represent less than 2 or more than 3 possible genotypes
	 * an empty array of size 3 is returned for the sample.
	 *
	 * If the BGEN probabilities represent 2 possible genotypes (biallelic for haploid samples)
	 * this is expanded as if diploid with a zero probability for heterozygosity.
	 *
	 * @param bgenProbabilities The BGEN probabilities returned by a BgenGenotypeData sampleVariantProvider.
	 * @return An array of arrays of size 3 with posterior probabilities.
	 */
	public static float[][] convertBgenProbabilitiesToProbabilities(double[][] bgenProbabilities) {
		// Define an array consisting of an array of posterior bgenProbabilities for each genotype
		float[][] probabilities = new float[bgenProbabilities.length][3];

		for (int sampleIndex = 0; sampleIndex < bgenProbabilities.length; sampleIndex++) {
			// Get the array of doubles
			double[] sampleProbabilitiesBgen = bgenProbabilities[sampleIndex];
			// Get the length of the probabilities array
			int probabilitiesArrayLength = sampleProbabilitiesBgen.length;
			// Initialize empty array for float values.
			float[] sampleProbabilities = new float[3];
			// Check for the length of the sample probabilities.
			if (probabilitiesArrayLength == 3) {
				// Convert the array of doubles to the array of floats
				IntStream.range(0, probabilitiesArrayLength)
						.forEach(index -> sampleProbabilities[index] = (float) sampleProbabilitiesBgen[index]);
			} else if (probabilitiesArrayLength == 2) {
				sampleProbabilities[0] = (float) sampleProbabilitiesBgen[0];
				sampleProbabilities[2] = (float) sampleProbabilitiesBgen[1];
			}
			// If probabilities array length is less than 2 or greater than 3,
			// just return an array of zeros. [0.0, 0.0, 0.0]

			// Insert the probabilities for this sample
			probabilities[sampleIndex] = sampleProbabilities;
		}

		return probabilities;
	}
}
