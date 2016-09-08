package org.molgenis.genotype.modifiable;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.util.FixedSizeIterable;
import org.molgenis.genotype.util.Ld;
import org.molgenis.genotype.util.LdCalculator;
import org.molgenis.genotype.util.LdCalculatorException;
import org.molgenis.genotype.util.MafCalculator;
import org.molgenis.genotype.util.MafResult;
import org.molgenis.genotype.variant.AbstractGeneticVariant;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.GeneticVariantMeta;
import org.molgenis.genotype.variant.GenotypeRecord;
import org.molgenis.genotype.variant.id.GeneticVariantId;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;

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
