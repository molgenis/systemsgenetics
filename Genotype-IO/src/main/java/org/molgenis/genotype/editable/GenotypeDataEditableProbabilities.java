/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.editable;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import org.apache.log4j.Logger;
import org.molgenis.genotype.AbstractRandomAccessGenotypeData;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.Sequence;
import org.molgenis.genotype.SimpleSequence;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.util.CalledDosageConvertor;
import org.molgenis.genotype.util.FixedSizeIterable;
import org.molgenis.genotype.util.ProbabilitiesConvertor;
import org.molgenis.genotype.util.RecordIteratorCreators;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.GeneticVariantMeta;
import org.molgenis.genotype.variant.GeneticVariantMetaMap;
import org.molgenis.genotype.variant.GenotypeRecord;
import org.molgenis.genotype.variant.ReadOnlyGeneticVariant;
import org.molgenis.genotype.variant.range.GeneticVariantRange;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantUniqueIdProvider;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;

/**
 *
 * @author Patrick Deelen
 */
public class GenotypeDataEditableProbabilities extends AbstractRandomAccessGenotypeData implements SampleVariantsProvider{

	/**
	 * variant, sample, probabilities
	 */
	private final float[][][] probabilitiesMatrix;
	
	private final int sampleVariantProviderUniqueId;
	private final GeneticVariantRange variants;
	private final HashMap<VariantInformation, Integer> variantAllelesIndex;
	private final LinkedHashMap<Sample, Integer> samples;
	private static final Logger LOGGER = Logger.getLogger(GenotypeDataEditableProbabilities.class);
	private final double minimumPosteriorProbabilityToCall;
	private final List<Boolean> phasing;
	private final LinkedHashSet<String> sequenceNames;
	private final GeneticVariantMeta geneticVariantMeta = GeneticVariantMetaMap.getGeneticVariantMetaGp();
	
	private static final double DEFAULT_MINIMUM_POSTERIOR_PROBABILITY_TO_CALL = 0.4f;

	public GenotypeDataEditableProbabilities(HashMap<VariantInformation, Integer> variantAllelesIndex, LinkedHashSet<Sample> samples) {
		this(variantAllelesIndex, samples, DEFAULT_MINIMUM_POSTERIOR_PROBABILITY_TO_CALL);
	}
	
	public GenotypeDataEditableProbabilities(HashMap<VariantInformation, Integer> variantAllelesIndex, LinkedHashSet<Sample> samples, double minimumPosteriorProbabilityToCall) {
		
		
		this.minimumPosteriorProbabilityToCall = minimumPosteriorProbabilityToCall;
		this.variantAllelesIndex = variantAllelesIndex;
		
		this.samples = new LinkedHashMap<Sample, Integer>(samples.size());
		int sampleCount = 0;
		for(Sample sample : samples){
			this.samples.put(sample, sampleCount);
			++sampleCount;
		}
			
		sampleVariantProviderUniqueId = SampleVariantUniqueIdProvider.getNextUniqueId();
		
		GeneticVariantRange.GeneticVariantRangeCreate variantFactory = GeneticVariantRange.createRangeFactory(variantAllelesIndex.size());
		
		this.sequenceNames = new LinkedHashSet<String>();
		for(VariantInformation variantInfo : variantAllelesIndex.keySet()){
			
			if(variantInfo.getAlleles().getAlleleCount() != 2){
				throw new GenotypeDataException("Only biallelic supported");
			}
			
			GeneticVariant variant = ReadOnlyGeneticVariant.createVariant(geneticVariantMeta, variantInfo.getVariantId(), variantInfo.getStartPos(), variantInfo.getSequenceName(), this, variantInfo.getAlleles().get(0), variantInfo.getAlleles().get(1));
			variantFactory.addVariant(variant);
			if(!sequenceNames.contains(variant.getSequenceName())){
				sequenceNames.add(variant.getSequenceName());
			}
			
		}
		
		this.variants = variantFactory.createRange();
		
		probabilitiesMatrix = new float[variants.size()][samples.size()][3];
		
		phasing = Collections.unmodifiableList(Collections.nCopies((int) samples.size(), false));
		
	}

	@Override
	public List<Sample> getSamples() {
		return new ArrayList<Sample>( samples.keySet() );
	}

	@Override
	public Map<String, Annotation> getVariantAnnotationsMap() {
		return Collections.emptyMap();
	}

	@Override
	public Map<String, SampleAnnotation> getSampleAnnotationsMap() {
		return Collections.emptyMap();
	}

	@Override
	public boolean isOnlyContaingSaveProbabilityGenotypes() {
		return true;
	}

	@Override
	public void close() throws IOException {
	}

	@Override
	public List<String> getSeqNames() {
		return new ArrayList<String>(sequenceNames);
	}

	@Override
	public Iterable<Sequence> getSequences() {
		List<Sequence> sequences = new ArrayList<Sequence>();
		for (String seqName : getSeqNames()) {
			sequences.add(new SimpleSequence(seqName, null, this));
		}

		return sequences;
	}

	@Override
	public Iterable<GeneticVariant> getVariantsByPos(String seqName, int startPos) {
		return variants.getVariantAtPos(seqName, startPos);
	}

	@Override
	public Iterable<GeneticVariant> getSequenceGeneticVariants(String seqName) {
		return variants.getVariantsBySequence(seqName);
	}

	@Override
	public Iterable<GeneticVariant> getVariantsByRange(String seqName, int rangeStart, int rangeEnd) {
		return variants.getVariantsByRange(seqName, rangeStart, rangeEnd);
	}

	@Override
	public List<Alleles> getSampleVariants(GeneticVariant variant) {
		return ProbabilitiesConvertor.convertProbabilitiesToAlleles(variant.getSampleGenotypeProbilities(), variant.getVariantAlleles(), minimumPosteriorProbabilityToCall);
	}

	@Override
	public List<Boolean> getSamplePhasing(GeneticVariant variant) {
		return phasing;
	}

	@Override
	public int cacheSize() {
		return variants.size();
	}

	@Override
	public int getSampleVariantProviderUniqueId() {
		return sampleVariantProviderUniqueId;
	}

	@Override
	public byte[] getSampleCalledDosage(GeneticVariant variant) {
		return CalledDosageConvertor.convertCalledAllelesToCalledDosage(variant.getSampleVariants(), variant.getVariantAlleles(), null);
	}

	@Override
	public float[] getSampleDosage(GeneticVariant variant) {
		return ProbabilitiesConvertor.convertProbabilitiesToDosage(variant.getSampleGenotypeProbilities(), minimumPosteriorProbabilityToCall);
	}

	@Override
	public float[][] getSampleProbilities(GeneticVariant variant) {
		return probabilitiesMatrix[variantAllelesIndex.get(new VariantInformation(variant.getPrimaryVariantId(), variant.getStartPos(), variant.getSequenceName(), variant.getVariantAlleles()))];
	}

	@Override
	public double[][] getSampleProbabilitiesComplex(GeneticVariant variant) {
		return ProbabilitiesConvertor.convertProbabilitiesToBgenProbabilities(getSampleProbilities(variant));
	}

	@Override
	public double[][][] getSampleProbabilitiesPhased(GeneticVariant variant) {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	public void setSampleVariantProbabilities(VariantInformation variantInfo, Sample sample, float[] probabilities){
		
		if(!samples.containsKey(sample)){
			throw new GenotypeDataException("Cannot set probabilities, sample not found: " + sample.getId() + " fam: " + sample.getFamilyId());
		}
		
		if(!variantAllelesIndex.containsKey(variantInfo)){
			throw new GenotypeDataException("Cannot set probabilities, variant not found: " + variantInfo.getSequenceName() + ":" + variantInfo.getStartPos());
		}
		
		if(probabilities.length != 3){
			throw new GenotypeDataException("Cannot set probabilities. Length not is 3");
		}
		
		probabilitiesMatrix[variantAllelesIndex.get(variantInfo)][samples.get(sample)] = probabilities;
		
	}
	
	@Override
	public FixedSizeIterable<GenotypeRecord> getSampleGenotypeRecords(GeneticVariant variant) {
		
		return RecordIteratorCreators.createIteratorFromProbs(variant.getSampleGenotypeProbilities());
		
	}
	
}
