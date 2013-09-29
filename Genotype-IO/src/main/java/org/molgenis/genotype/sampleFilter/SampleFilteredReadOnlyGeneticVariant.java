/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.sampleFilter;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.probabilities.SampleVariantProbabilities;
import org.molgenis.genotype.util.Ld;
import org.molgenis.genotype.util.LdCalculatorException;
import org.molgenis.genotype.variant.AbstractGeneticVariant;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.id.GeneticVariantId;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;

/**
 *
 * @author Patrick Deelen
 */
public class SampleFilteredReadOnlyGeneticVariant extends AbstractGeneticVariant {

	private final GeneticVariant original;
	private final SampleFilterableGenotypeData genotypeData;

	public SampleFilteredReadOnlyGeneticVariant(GeneticVariant original, SampleFilterableGenotypeData genotypeData) {
		this.original = original;
		this.genotypeData = genotypeData;
	}

	@Override
	public String getPrimaryVariantId() {
		return original.getPrimaryVariantId();
	}

	@Override
	public List<String> getAlternativeVariantIds() {
		return original.getAlternativeVariantIds();
	}

	@Override
	public List<String> getAllIds() {
		return original.getAllIds();
	}

	@Override
	public GeneticVariantId getVariantId() {
		return original.getVariantId();
	}

	@Override
	public int getStartPos() {
		return original.getStartPos();
	}

	@Override
	public String getSequenceName() {
		return original.getSequenceName();
	}

	@Override
	public Alleles getVariantAlleles() {
		return original.getVariantAlleles();
	}

	@Override
	public int getAlleleCount() {
		return original.getAlleleCount();
	}

	@Override
	public Allele getRefAllele() {
		return original.getRefAllele();
	}

	@Override
	public List<Alleles> getSampleVariants() {


		ArrayList<Alleles> sampleAlleles = new ArrayList<Alleles>(original.getSampleVariants());

		Iterator<Alleles> sampleAllelesIterator = sampleAlleles.iterator();
		Iterator<Sample> sampleIterator = genotypeData.getOriginalSampleList().iterator();

		try {
			while (sampleAllelesIterator.hasNext()) {
				sampleAllelesIterator.next();
				if (!genotypeData.getSampleFilter().doesSamplePassFilter(sampleIterator.next())) {
					sampleAllelesIterator.remove();
				}
			}
		} catch (NoSuchElementException e) {
			throw new GenotypeDataException("Error in filtering on included samples. More alleles than samples detected", e);
		}

		return Collections.unmodifiableList(sampleAlleles);

	}

	@Override
	public Map<String, ?> getAnnotationValues() {
		return original.getAnnotationValues();
	}

	@Override
	public double getMinorAlleleFrequency() {
		return original.getMinorAlleleFrequency();
	}

	@Override
	public Allele getMinorAllele() {
		return original.getMinorAllele();
	}

	@Override
	public boolean isSnp() {
		return original.isSnp();
	}

	@Override
	public boolean isAtOrGcSnp() {
		return original.isAtOrGcSnp();
	}

	@Override
	public Ld calculateLd(GeneticVariant other) throws LdCalculatorException {
		return original.calculateLd(other);
	}

	@Override
	public boolean isBiallelic() {
		return original.isBiallelic();
	}

	@Override
	public float[] getSampleDosages() {

		float[] unfilteredDosages = original.getSampleDosages();
		float[] includedSamplesDosages = new float[genotypeData.getIncludedSampleCount()];

		Iterator<Sample> sampleIterator = genotypeData.getOriginalSampleList().iterator();

		try {
			int i = 0;
			for (float dosage : unfilteredDosages) {
				if (genotypeData.getSampleFilter().doesSamplePassFilter(sampleIterator.next())) {
					includedSamplesDosages[i] = dosage;
					++i;
				}
			}
		} catch (NoSuchElementException e) {
			throw new GenotypeDataException("Error in filtering on included samples. More dosage values than samples detected", e);
		}

		return includedSamplesDosages;

	}

	@Override
	public byte[] getSampleCalledDosages() {

		byte[] unfilteredDosages = original.getSampleCalledDosages();
		byte[] includedSamplesDosages = new byte[genotypeData.getIncludedSampleCount()];

		Iterator<Sample> sampleIterator = genotypeData.getOriginalSampleList().iterator();

		try {
			int i = 0;
			for (byte dosage : unfilteredDosages) {
				if (genotypeData.getSampleFilter().doesSamplePassFilter(sampleIterator.next())) {
					includedSamplesDosages[i] = dosage;
					++i;
				}
			}
		} catch (NoSuchElementException e) {
			throw new GenotypeDataException("Error in filtering on included samples. More dosage values than samples detected", e);
		}

		return includedSamplesDosages;

	}

	@Override
	public List<Boolean> getSamplePhasing() {

		ArrayList<Boolean> samplePhasing = new ArrayList<Boolean>(original.getSamplePhasing());

		Iterator<Boolean> samplePhasingIterator = samplePhasing.iterator();
		Iterator<Sample> sampleIterator = genotypeData.getOriginalSampleList().iterator();

		try {
			while (samplePhasingIterator.hasNext()) {
				samplePhasingIterator.next();
				if (!genotypeData.getSampleFilter().doesSamplePassFilter(sampleIterator.next())) {
					samplePhasingIterator.remove();
				}
			}
		} catch (NoSuchElementException e) {
			throw new GenotypeDataException("Error in filtering on included samples. More alleles than samples detected", e);
		}

		return Collections.unmodifiableList(samplePhasing);

	}

	@Override
	public SampleVariantProbabilities[] getSampleGenotypeProbilities() {
		
		SampleVariantProbabilities[] unfilteredProbs = original.getSampleGenotypeProbilities();
		SampleVariantProbabilities[] includedSamplesProbs = new SampleVariantProbabilities[genotypeData.getIncludedSampleCount()];

		Iterator<Sample> sampleIterator = genotypeData.getOriginalSampleList().iterator();

		try {
			int i = 0;
			for (SampleVariantProbabilities prob : unfilteredProbs) {
				if (genotypeData.getSampleFilter().doesSamplePassFilter(sampleIterator.next())) {
					includedSamplesProbs[i] = prob;
					++i;
				}
			}
		} catch (NoSuchElementException e) {
			throw new GenotypeDataException("Error in filtering on included samples. More prob values than samples detected", e);
		}

		return includedSamplesProbs;
		
	}
	
	

	@Override
	public SampleVariantsProvider getSampleVariantsProvider() {
		return original.getSampleVariantsProvider();
	}
}
