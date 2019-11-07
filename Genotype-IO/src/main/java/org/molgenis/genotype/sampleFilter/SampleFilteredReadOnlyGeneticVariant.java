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
import org.molgenis.genotype.util.*;
import org.molgenis.genotype.variant.AbstractGeneticVariant;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.GeneticVariantMeta;
import org.molgenis.genotype.variant.GenotypeRecord;
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
	public Alleles getAlternativeAlleles() {
		return original.getAlternativeAlleles();
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
		// Do not cache MAF results since sample filter might change
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
		// Do not cache MAF results since sample filter might change
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
	public float[][] getSampleGenotypeProbilities() {

		float[][] unfilteredProbs = original.getSampleGenotypeProbilities();
		float[][] includedSamplesProbs = new float[genotypeData.getIncludedSampleCount()][3];

		Iterator<Sample> sampleIterator = genotypeData.getOriginalSampleList().iterator();

		try {
			int i = 0;
			for (float[] prob : unfilteredProbs) {
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
	public double[][] getSampleGenotypeProbabilitiesBgen() {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	@Override
	public SampleVariantsProvider getSampleVariantsProvider() {
		return original.getSampleVariantsProvider();
	}

	@Override
	public GeneticVariantMeta getVariantMeta() {
		return original.getVariantMeta();
	}

	@Override
	public FixedSizeIterable<GenotypeRecord> getSampleGenotypeRecords() {

		final Iterator<GenotypeRecord> originalRecords = original.getSampleGenotypeRecords().iterator();

		return new FixedSizeIterable<GenotypeRecord>() {
			@Override
			public int size() {
				return genotypeData.getIncludedSampleCount();
			}

			@Override
			public Iterator<GenotypeRecord> iterator() {
				return new Iterator<GenotypeRecord>() {
					int i = 0;
					Iterator<Sample> sampleIterator = genotypeData.getOriginalSampleList().iterator();

					@Override
					public boolean hasNext() {
						return i < genotypeData.getIncludedSampleCount();
					}

					@Override
					public GenotypeRecord next() {
						if (i < genotypeData.getIncludedSampleCount() ) {
							try {
								while (originalRecords.hasNext()) {
									GenotypeRecord x = originalRecords.next();
									if (genotypeData.getSampleFilter().doesSamplePassFilter(sampleIterator.next())) {
										++i;
										return x;
									}
								}
								throw new GenotypeDataException("Error in filtering on included samples.");
							} catch (NoSuchElementException e) {
								throw new GenotypeDataException("Error in filtering on included samples. More records than samples detected", e);
							}
						} else {
							throw new GenotypeDataException("Error in filtering on included samples.");
						}
					}

					@Override
					public void remove() {
						throw new UnsupportedOperationException("Not supported yet.");
					}
				};
			}
		};


	}
}