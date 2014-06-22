/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.editable;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import org.apache.log4j.Logger;
import org.molgenis.genotype.AbstractRandomAccessGenotypeData;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.Sequence;
import org.molgenis.genotype.SimpleSequence;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.ReadOnlyGeneticVariant;
import org.molgenis.genotype.variant.range.GeneticVariantRange;

/**
 *
 * @author Patrick Deelen
 */
public class GenotypeDataCustomVariantProvider<S extends EditableSampleVariantsProvider> extends AbstractRandomAccessGenotypeData {

	private final GeneticVariantRange variants;
	private final List samples;
	private static final Logger LOGGER = Logger.getLogger(GenotypeDataCustomVariantProvider.class);	
	private final LinkedHashSet<String> sequenceNames;
	private final S sampleVariantProvider;

	public GenotypeDataCustomVariantProvider(Iterable<GeneticVariant> variants, List<Sample> samples, S sampleVariantProvider) {

		this.sampleVariantProvider = sampleVariantProvider;

		this.samples = samples;

		GeneticVariantRange.GeneticVariantRangeCreate variantFactory = GeneticVariantRange.createRangeFactory();

		this.sequenceNames = new LinkedHashSet<String>();
		for (GeneticVariant originalVariant : variants) {

			if (originalVariant.getVariantAlleles().getAlleleCount() != 2) {
				throw new GenotypeDataException("Only biallelic supported");
			}

			GeneticVariant variant = ReadOnlyGeneticVariant.createVariant(this.sampleVariantProvider.getGeneticVariantMeta(), originalVariant.getVariantId(), originalVariant.getStartPos(), originalVariant.getSequenceName(), this.sampleVariantProvider, originalVariant.getVariantAlleles().get(0), originalVariant.getVariantAlleles().get(1));
			variantFactory.addVariant(variant);
			if (!sequenceNames.contains(variant.getSequenceName())) {
				sequenceNames.add(variant.getSequenceName());
			}

		}

		this.variants = variantFactory.createRange();

		

	}

	@Override
	public List<Sample> getSamples() {
		return samples;
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

	public S getSampleVariantProvider() {
		return sampleVariantProvider;
	}

}
