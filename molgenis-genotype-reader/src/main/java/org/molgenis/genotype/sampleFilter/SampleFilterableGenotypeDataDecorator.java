/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.sampleFilter;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import org.molgenis.genotype.AbstractRandomAccessGenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.Sequence;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 *
 * @author Patrick Deelen
 */
public class SampleFilterableGenotypeDataDecorator extends AbstractRandomAccessGenotypeData implements SampleFilterableGenotypeData{
	
	private final RandomAccessGenotypeData original;
	private final ArrayList<Sample> includedSamples;
	private final SampleFilter sampleFilter;

	public SampleFilterableGenotypeDataDecorator(RandomAccessGenotypeData original, SampleFilter sampleFilter) {
		this.original = original;
		this.sampleFilter = sampleFilter;

		includedSamples = new ArrayList<Sample>();

		for (Sample sample : original.getSamples()) {
			if (sampleFilter.doesSamplePassFilter(sample)) {
				includedSamples.add(sample);
			}
		}

	}

	@Override
	public List<String> getSeqNames() {
		return original.getSeqNames();
	}

	@Override
	public Iterable<Sequence> getSequences() {
		return original.getSequences();
	}

	@Override
	public Sequence getSequenceByName(String name) {
		return original.getSequenceByName(name);
	}

	@Override
	public Iterable<GeneticVariant> getVariantsByPos(String seqName, int startPos) {
		return new ConvertToSampleFilterIterable(original.getVariantsByPos(seqName, startPos), this);
	}

	@Override
	public GeneticVariant getSnpVariantByPos(String seqName, int startPos) {
		GeneticVariant variant = original.getSnpVariantByPos(seqName, startPos);
		return variant == null ? null : new SampleFilteredReadOnlyGeneticVariant(variant, this);
	}

	@Override
	public Iterable<GeneticVariant> getSequenceGeneticVariants(String seqName) {
		return new ConvertToSampleFilterIterable(original.getSequenceGeneticVariants(seqName), this);
	}

	@Override
	public Iterable<GeneticVariant> getVariantsByRange(String seqName, int rangeStart, int rangeEnd) {
		return new ConvertToSampleFilterIterable(original.getVariantsByRange(seqName, rangeStart, rangeEnd), this);
	}

	@Override
	public List<Annotation> getVariantAnnotations() {
		return original.getVariantAnnotations();
	}

	@Override
	public Annotation getVariantAnnotation(String annotationId) {
		return original.getVariantAnnotation(annotationId);
	}

	@Override
	public List<SampleAnnotation> getSampleAnnotations() {
		return original.getSampleAnnotations();
	}

	@Override
	public Annotation getSampleAnnotation(String annotationId) {
		return original.getSampleAnnotation(annotationId);
	}

	@Override
	public List<Sample> getSamples() {
		return includedSamples;
	}

	@Override
	public Iterator<GeneticVariant> iterator() {
		return new ConvertToSampleFilterIterable(original, this).iterator();
	}

	@Override
	public int getIncludedSampleCount() {
		return includedSamples.size();
	}

	@Override
	public List<Sample> getOriginalSampleList() {
		return original.getSamples();
	}

	@Override
	public SampleFilter getSampleFilter() {
		return sampleFilter;
	}

	@Override
	public void close() throws IOException {
		original.close();
	}

	@Override
	public Map<String, SampleAnnotation> getSampleAnnotationsMap() {
		return original.getSampleAnnotationsMap();
	}

	@Override
	public Map<String, ? extends Annotation> getVariantAnnotationsMap() {
		return original.getVariantAnnotationsMap();
	}
	
}
