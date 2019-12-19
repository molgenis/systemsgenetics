package org.molgenis.genotype.modifiable;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.util.*;
import org.molgenis.genotype.variant.AbstractGeneticVariant;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.GeneticVariantMeta;
import org.molgenis.genotype.variant.GenotypeRecord;
import org.molgenis.genotype.variant.id.GeneticVariantId;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;

import static org.molgenis.genotype.bgen.BgenGenotypeWriter.getPloidy;
import static org.molgenis.genotype.bgen.BgenGenotypeData.getAlleleCountsPerProbability;

public class ModifiableGeneticVariant extends AbstractGeneticVariant {

	private final GeneticVariant originalVariant;
	private final ModifiableGenotypeData modifiableGenotypeData;

	public ModifiableGeneticVariant(GeneticVariant originalVariant, ModifiableGenotypeData modifiableGenotypeData) {
		this.originalVariant = originalVariant;
		this.modifiableGenotypeData = modifiableGenotypeData;
	}

	@Override
	public String getPrimaryVariantId() {
		return getVariantId().getPrimairyId();
	}

	@Override
	public List<String> getAlternativeVariantIds() {
		return getVariantId().getAlternativeIds();
	}

	@Override
	public List<String> getAllIds() {
		return getVariantId().getVariantIds();
	}

	@Override
	public GeneticVariantId getVariantId() {
		GeneticVariantId variantId = modifiableGenotypeData.getUpdatedId(this);
		if (variantId != null) {
			return variantId;
		} else {
			return originalVariant.getVariantId();
		}
	}

	@Override
	public int getStartPos() {
		return originalVariant.getStartPos();
	}

	@Override
	public String getSequenceName() {
		return originalVariant.getSequenceName();
	}

	@Override
	public Alleles getVariantAlleles() {
		Alleles alleles = modifiableGenotypeData.getUpdatedAlleles(this);
		if (alleles != null) {
			return alleles;
		} else {
			return originalVariant.getVariantAlleles();
		}
	}

	@Override
	public Alleles getAlternativeAlleles() {
		ArrayList<Allele> altAlleles = new ArrayList<>(this.getVariantAlleles().getAlleles());
		altAlleles.remove(this.getRefAllele());
		return Alleles.createAlleles(altAlleles);
	}

	@Override
	public int getAlleleCount() {
		return getVariantAlleles().getAlleleCount();
	}

	@Override
	public Allele getRefAllele() {
		Allele refAllele = modifiableGenotypeData.getUpdatedRef(this);
		if (refAllele != null) {
			return refAllele;
		} else {
			return originalVariant.getRefAllele();
		}
	}

	@Override
	public List<Alleles> getSampleVariants() {
		return getSampleVariantsProvider().getSampleVariants(originalVariant);
	}

	@Override
	public List<Boolean> getSamplePhasing() {
		return getSampleVariantsProvider().getSamplePhasing(originalVariant);
	}

	@Override
	public Map<String, ?> getAnnotationValues() {
		return originalVariant.getAnnotationValues();
	}

	@Override
	public double getMinorAlleleFrequency() {
		// Do not cache MAF results since modifications to alleles need to be
		// reflected
		try {
			MafResult mafResult = MafCalculator.calculateMaf(getVariantAlleles(), getRefAllele(), getSampleVariants());
			return mafResult.getFreq();
		} catch (NullPointerException e) {
			throw new GenotypeDataException("NullPointerException in maf caculation. " + getVariantAlleles() + " ref: "
					+ getRefAllele(), e);
		}

	}

	@Override
	public Allele getMinorAllele() {
		// Do not cache MAF results since modifications to alleles need to be
		// reflected
		try {
			MafResult mafResult = MafCalculator.calculateMaf(getVariantAlleles(), getRefAllele(), getSampleVariants());
			return mafResult.getMinorAllele();
		} catch (NullPointerException e) {
			throw new GenotypeDataException("NullPointerException in maf caculation. " + getVariantAlleles() + " ref: "
					+ getRefAllele(), e);
		}
	}

	@Override
	public boolean isSnp() {
		return getVariantAlleles().isSnp();
	}

	@Override
	public boolean isAtOrGcSnp() {
		return getVariantAlleles().isAtOrGcSnp();
	}

	@Override
	public Ld calculateLd(GeneticVariant other) throws LdCalculatorException {
		return LdCalculator.calculateLd(this, other);
	}

	@Override
	public boolean isBiallelic() {
		return getVariantAlleles().getAlleleCount() == 2;
	}

	@Override
	public float[] getSampleDosages() {
		float[] dosageByProvider = getSampleVariantsProvider().getSampleDosage(originalVariant);

		Allele refUsedForOriginalDosage = originalVariant.getRefAllele() == null ? originalVariant.getVariantAlleles()
				.get(0) : originalVariant.getRefAllele();

		if(modifiableGenotypeData.isSwapped(originalVariant)){
			refUsedForOriginalDosage = refUsedForOriginalDosage.getComplement();
		}
		
		Allele refShouldBeUsed = getRefAllele() == null ? getVariantAlleles().get(0) : getRefAllele();
		
		
		if (refUsedForOriginalDosage == refShouldBeUsed) {
			return dosageByProvider;
		}  else {

			// Here we have to do the swap of the dosage to match the new ref.
			// (org - 2 ) * 1
			float[] newDosage = new float[dosageByProvider.length];
			for (int i = 0; i < dosageByProvider.length; ++i) {
				// -1 -> 1 -> -1
				// 0 -> 0 -> 2
				// 0.5 -> -0.5 -> 1.5
				// 1 -> -1 -> 1
				// 1.5 -> -1.5 -> 0.5
				// 2 -> -2 -> 0
				newDosage[i] = dosageByProvider[i] == -1 ? -1 : (dosageByProvider[i] * -1) + 2;
			}
			return newDosage;
		}

	}

	@Override
	public byte[] getSampleCalledDosages() {

		byte[] dosageByProvider = getSampleVariantsProvider().getSampleCalledDosage(originalVariant);

		Allele refUsedForOriginalDosage = originalVariant.getRefAllele() == null ? originalVariant.getVariantAlleles()
				.get(0) : originalVariant.getRefAllele();

		if(modifiableGenotypeData.isSwapped(originalVariant)){
			refUsedForOriginalDosage = refUsedForOriginalDosage.getComplement();
		}
		

		Allele refShouldBeUsed = getRefAllele() == null ? getVariantAlleles().get(0) : getRefAllele();

		// System.out.println("Ref used: " + refUsedForOriginalDosage +
		// " should be: " + refShouldBeUsed);

		if (refUsedForOriginalDosage == refShouldBeUsed) {
			return dosageByProvider;
		} else {

			// Here we have to do the swap of the dosage to match the new ref.
			// (org - 2 ) * 1
			byte[] newDosage = new byte[dosageByProvider.length];
			for (int i = 0; i < dosageByProvider.length; ++i) {
				// -1 -> 1 -> -1
				// 0 -> 0 -> 2
				// 0.5 -> -0.5 -> 1.5
				// 1 -> -1 -> 1
				// 1.5 -> -1.5 -> 0.5
				// 2 -> -2 -> 0
				newDosage[i] = (byte) (dosageByProvider[i] == -1 ? -1 : (dosageByProvider[i] * -1) + 2);
			}
			return newDosage;
		}

	}

	@Override
	public float[][] getSampleGenotypeProbilities() {

		float[][] probByProvider = getSampleVariantsProvider().getSampleProbilities(originalVariant);

		Allele refUsedForOriginalDosage = originalVariant.getRefAllele() == null ? originalVariant.getVariantAlleles()
				.get(0) : originalVariant.getRefAllele();

		if(modifiableGenotypeData.isSwapped(originalVariant)){
			refUsedForOriginalDosage = refUsedForOriginalDosage.getComplement();
		}
		

		Allele refShouldBeUsed = getRefAllele() == null ? getVariantAlleles().get(0) : getRefAllele();
		
		
		if (refUsedForOriginalDosage == refShouldBeUsed) {
			return probByProvider;
		} else {

			float[][] probs = new float[probByProvider.length][3];

			for (int i = 0; i < probByProvider.length; ++i) {

				for (int j = 0; j < 3; ++j) {
					probs[i][2 - j] = probByProvider[i][j];
				}

			}

			return probs;

		}

	}

	@Override
	public double[][] getSampleGenotypeProbabilitiesComplex() {

		double[][] probByProvider = getSampleVariantsProvider().getSampleProbabilitiesComplex(originalVariant);

		Allele refUsedForOriginalDosage = originalVariant.getRefAllele() == null ? originalVariant.getVariantAlleles()
				.get(0) : originalVariant.getRefAllele();

		if(modifiableGenotypeData.isSwapped(originalVariant)){
			refUsedForOriginalDosage = refUsedForOriginalDosage.getComplement();
		}

		Allele refShouldBeUsed = getRefAllele() == null ? getVariantAlleles().get(0) : getRefAllele();

		if (refUsedForOriginalDosage == refShouldBeUsed) {
			return probByProvider;
		} else {
			// Currently throw an exception whenever the following is not true
			if (!(isBiallelic() && Arrays.stream(probByProvider).allMatch(a -> a.length == 3))) {
				throw new GenotypeDataException("Can't swap probabilities multiallelic / polyploid variants");
			}
			// Here we have to swap the probabilities to match the new ref
			// Probabilities for complex probabilities are sorted in a special manner
			// We can get a list of allele counts corresponding to the genotypes in current order,
			// which can be easily ordered into the situation with the new ref allele.
			// Therefore we can sort the probabilities using the list of allele counts as well.

			// First generate a list with alleles, in which every allele is encoded as the index in the new situation
			// Get the index of the new reference allele within the old list of alleles
			// TODO: implement obtaining the index of the new reference allele according to the old list of alleles.
			int refShouldBeUsedOldIndex = 1; // first index == second allele
			// Continue to create a list of alleles
			List<Integer> allelesAsIntegers = IntStream.range(1, getAlleleCount()).boxed().collect(Collectors.toList());
			allelesAsIntegers.add(refShouldBeUsedOldIndex, 0); // Add the ref with index 0;

			double[][] probs = new double[probByProvider.length][];

			for (int i = 0; i < probByProvider.length; ++i) {
				// Get the number of probabilities for this sample
				double[] sampleProbabilities = probByProvider[i];
				int numberOfProbabilities = sampleProbabilities.length;

				// Convert the array to a list to be able to make use of the .indexOf(...) method.
				List<Pair<Integer, Double>> probabilitiesWithIndices = new ArrayList<>();
				double[] doubles = probByProvider[i];
				for (int j = 0; j < doubles.length; j++) {
					double probability = doubles[j];
					probabilitiesWithIndices.add(new ImmutablePair<>(j, probability));
				}
				// Get the allele counts per probability, in which the alleles are sorted according to the old order alleles.
				List<List<Integer>> alleleCounts = getAlleleCountsPerProbability(
						allelesAsIntegers,
						getPloidy(numberOfProbabilities, getAlleleCount()));
				// Now sort the probabilities for this sample
				probabilitiesWithIndices.sort((left, right) -> {
					// Initialize return value at 0 (no difference)
					int comparison = 0;
					// Loop through the allele counts from the last to first to match sorting precedence.
					// The first column (from last to first) that shows difference indicates sorting order,
					// so continue until the comparison value is not equal to zero.
					for (int alleleIndex = alleleCounts.get(left.getLeft()).size() - 1; // Initialize at last index
						 comparison == 0 && alleleIndex >= 0; // Stay in loop condition
						 alleleIndex--) { // Decrease index
						comparison = alleleCounts.get(left.getLeft()).get(alleleIndex)
								.compareTo(
										alleleCounts.get(right.getLeft()).get(alleleIndex));
					}
					return comparison;
				});
				// Insert the sorted probabilities as an array into the new probabilities.
				probs[i] = probabilitiesWithIndices.stream().mapToDouble(Pair::getRight).toArray();

			}
			return probs;

		}
	}

	@Override
	public double[][][] getSampleGenotypeProbabilitiesPhased() {
		double[][][] probByProvider = getSampleVariantsProvider().getSampleProbabilitiesPhased(originalVariant);

		Allele refUsedForOriginalDosage = originalVariant.getRefAllele() == null ? originalVariant.getVariantAlleles()
				.get(0) : originalVariant.getRefAllele();

		if(modifiableGenotypeData.isSwapped(originalVariant)){
			refUsedForOriginalDosage = refUsedForOriginalDosage.getComplement();
		}

		Allele refShouldBeUsed = getRefAllele() == null ? getVariantAlleles().get(0) : getRefAllele();

		if (refUsedForOriginalDosage == refShouldBeUsed) {
			return probByProvider;
		} else {
			// Currently throw an exception whenever the following is not true
			if (!(isBiallelic() && Arrays.stream(probByProvider).allMatch(a -> a.length == 2))) {
				throw new GenotypeDataException("Can't swap probabilities multiallelic / polyploid variants");
			}

			double[][][] probs = new double[probByProvider.length][2][2];

			for (int i = 0; i < probByProvider.length; i++) {
				for (int j = 0; j < 2; j++) {
					probs[i][j][1] = probByProvider[i][j][0];
					probs[i][j][0] = probByProvider[i][j][1];
				}
			}
			return probs;
		}
	}

	@Override
	public SampleVariantsProvider getSampleVariantsProvider() {
		SampleVariantsProvider sampleVariantProvider = modifiableGenotypeData.getUpdatedSampleVariantProvider(this);
		if (sampleVariantProvider != null) {
			return sampleVariantProvider;
		} else {
			return originalVariant.getSampleVariantsProvider();
		}
	}

	/**
	 * @return the originalVariant
	 */
	protected GeneticVariant getOriginalVariant() {
		return originalVariant;
	}

	/**
	 * Updates reference allele
	 *
	 * @param newRefAllele
	 * @throws GenotypeModificationException if reference allele is not one of
	 * the variant alleles
	 */
	public void updateRefAllele(Allele newRefAllele) {
		modifiableGenotypeData.updateRefAllele(this, newRefAllele);
	}

	/**
	 * Updates reference allele
	 *
	 * @param newRefAllele
	 * @throws GenotypeModificationException if reference allele is not one of
	 * the variant alleles
	 */
	public void updateRefAllele(String newRefAllele) {
		updateRefAllele(Allele.create(newRefAllele));
	}

	/**
	 * Updates reference allele
	 *
	 * @param newRefAllele
	 * @throws GenotypeModificationException if reference allele is not one of
	 * the variant alleles
	 */
	public void updateRefAllele(char newRefAllele) {
		updateRefAllele(Allele.create(newRefAllele));
	}

	/**
	 * Sets new primary ID. Old ID will made an alternative ID
	 *
	 * @param newPrimaryId
	 */
	public void updatePrimaryId(String newPrimaryId) {
		modifiableGenotypeData.updateVariantPrimaryId(this, newPrimaryId);
	}

	/**
	 * Overwrite old variant ID
	 *
	 * @param newVariantId
	 */
	public void updateId(GeneticVariantId newVariantId) {
		modifiableGenotypeData.updateVariantId(this, newVariantId);
	}

	/**
	 * Swap the alleles from this variants. Variant Alleles, Reference Allele
	 * and Sample Alleles will be updated
	 */
	public void swap() {
		modifiableGenotypeData.swapGeneticVariant(this);
	}

	/**
	 * Exclude this variant from the modifiable genotype data it belongs to.
	 */
	public void exclude() {
		modifiableGenotypeData.excludeVariant(this);
	}

	@Override
	public GeneticVariantMeta getVariantMeta() {
		return originalVariant.getVariantMeta();
	}

	@Override
	public FixedSizeIterable<GenotypeRecord> getSampleGenotypeRecords() {
		throw new UnsupportedOperationException("Genotype records are currently not implemented for modifiable genetic variants.");
	}
}
