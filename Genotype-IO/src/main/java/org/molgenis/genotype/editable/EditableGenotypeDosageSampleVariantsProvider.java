package org.molgenis.genotype.editable;

import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.util.CalledDosageConvertor;
import org.molgenis.genotype.util.FixedSizeIterable;
import org.molgenis.genotype.util.ProbabilitiesConvertor;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.GeneticVariantMeta;
import org.molgenis.genotype.variant.GeneticVariantMetaMap;
import org.molgenis.genotype.variant.GenotypeRecord;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantUniqueIdProvider;

/**
 *
 * @author Patrick Deelen
 */
public class EditableGenotypeDosageSampleVariantsProvider implements EditableSampleVariantsProvider {
	
	private final int sampleVariantProviderUniqueId;
	private final GeneticVariantMeta geneticVariantMeta;
	private final int numberSamples;
	private final List<Boolean> phasing;
	private final HashMap<VariantInformation, List<Alleles>> sampleVariants;
	private final HashMap<VariantInformation, float[]> sampleDosages;

	public EditableGenotypeDosageSampleVariantsProvider(int numberSamples) {
		sampleVariantProviderUniqueId = SampleVariantUniqueIdProvider.getNextUniqueId();
		
		HashMap<String, GeneticVariantMeta.Type> variantMeta = new HashMap<String, GeneticVariantMeta.Type>(2);
		variantMeta.put("GT", GeneticVariantMeta.Type.ALLELES);
		variantMeta.put("DS", GeneticVariantMeta.Type.FLOAT);
		this.geneticVariantMeta = GeneticVariantMetaMap.createGeneticVariantMeta(variantMeta);
		
		this.numberSamples = numberSamples;
		phasing = Collections.unmodifiableList(Collections.nCopies(numberSamples, false));
		
		sampleVariants = new HashMap<VariantInformation, List<Alleles>>();
		sampleDosages = new HashMap<VariantInformation, float[]>();
		
	}
	

	@Override
	public List<Alleles> getSampleVariants(GeneticVariant variant) {
		List<Alleles> alleles = sampleVariants.get(new VariantInformation(variant.getPrimaryVariantId(), variant.getStartPos(), variant.getSequenceName(), variant.getVariantAlleles()));
		if(alleles == null){
			throw new GenotypeDataException("No alleles for variant: " + variant.getPrimaryVariantId() + " in editable data");
		}
		return alleles;
	}

	@Override
	public FixedSizeIterable<GenotypeRecord> getSampleGenotypeRecords(GeneticVariant variant) {
		return new FixedSizeIterableTriTyper(variant);
	}

	@Override
	public List<Boolean> getSamplePhasing(GeneticVariant variant) {
		return phasing;
	}

	@Override
	public boolean arePhasedProbabilitiesPresent(GeneticVariant variant) {
		return false;
	}

	@Override
	public int cacheSize() {
		return sampleVariants.size();
	}

	@Override
	public int getSampleVariantProviderUniqueId() {
		return sampleVariantProviderUniqueId;
	}

	@Override
	public byte[] getSampleCalledDosage(GeneticVariant variant) {
		return CalledDosageConvertor.convertCalledAllelesToCalledDosage(this.getSampleVariants(variant),
				variant.getVariantAlleles(), variant.getRefAllele());
	}

	@Override
	public float[] getSampleDosage(GeneticVariant variant) {
		float[] dosages = sampleDosages.get(new VariantInformation(variant.getPrimaryVariantId(), variant.getStartPos(), variant.getSequenceName(), variant.getVariantAlleles()));
		if(dosages == null){
			throw new GenotypeDataException("No dosage for variant: " + variant.getPrimaryVariantId() + " in editable data");
		}
		return dosages;
	}

	@Override
	public float[][] getSampleProbilities(GeneticVariant variant) {
		return ProbabilitiesConvertor.convertDosageToProbabilityHeuristic(variant.getSampleDosages());
	}

	@Override
	public double[][] getSampleProbabilitiesComplex(GeneticVariant variant) {
		return ProbabilitiesConvertor.convertProbabilitiesToComplexProbabilities(getSampleProbilities(variant));
	}

	@Override
	public double[][][] getSampleProbabilitiesPhased(GeneticVariant variant) {
		throw new GenotypeDataException("Phased data not available");
	}

	@Override
	public GeneticVariantMeta getGeneticVariantMeta() {
		return this.geneticVariantMeta;
	}
	
	public void setDosageAndGenotypes(VariantInformation variant, List<Alleles> alleles, float[] dosages){
		if(alleles.size() != dosages.length){
			throw new GenotypeDataException();
		}
		
		if(alleles.size() != numberSamples){
			throw new GenotypeDataException();
		}
		
		sampleVariants.put(variant, alleles);
		sampleDosages.put(variant, dosages);
		
	}
	
	public void setDosageAndGenotypes(GeneticVariant variant, List<Alleles> alleles, float[] dosages){
		this.setDosageAndGenotypes(new VariantInformation(variant.getPrimaryVariantId(), variant.getStartPos(), variant.getSequenceName(), variant.getVariantAlleles()), alleles, dosages);
	}
	
	private class FixedSizeIterableTriTyper implements FixedSizeIterable<GenotypeRecord> {

		private final List<Alleles> alleles;
		private final float[] dosages;

		public FixedSizeIterableTriTyper(GeneticVariant variant) {
			alleles = variant.getSampleVariants();
			dosages = variant.getSampleDosages();
		}

		@Override
		public int size() {
			return alleles.size();
		}

		@Override
		public Iterator<GenotypeRecord> iterator() {
			return new Iterator<GenotypeRecord>() {
				int i = 0;

				@Override
				public boolean hasNext() {
					return i + 1 < alleles.size();
				}

				@Override
				public GenotypeRecord next() {
					++i;
					return new FixedSizeIterableTriTyper.TriTyperGenotypeRecord(i);
				}

				@Override
				public void remove() {
					throw new UnsupportedOperationException("Not supported ever yet.");
				}
			};
		}

		private class TriTyperGenotypeRecord implements GenotypeRecord {

			private final int i;

			public TriTyperGenotypeRecord(int i) {
				this.i = i;
			}

			@Override
			public Object getGenotypeRecordData(String recordId) {
				if (recordId.equals("GT")) {
					return alleles.get(i);
				} else if (recordId.equals("DS")) {
					return dosages[i];
				} else {
					return null;
				}
			}

			@Override
			public Alleles getSampleAlleles() {
				return alleles.get(i);
			}

			@Override
			public float[] getSampleProbs() {
				return null;
			}

			@Override
			public float getSampleDosage() {
				return dosages[i];
			}

			@Override
			public boolean containsGenotypeRecord(String recordId) {
				return recordId.equals("GT") | recordId.equals("DS");
			}
		}
	}

}
