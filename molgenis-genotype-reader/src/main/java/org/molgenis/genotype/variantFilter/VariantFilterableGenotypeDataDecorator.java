/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.variantFilter;

import java.io.IOException;
import java.util.Iterator;
import java.util.List;
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
public class VariantFilterableGenotypeDataDecorator implements RandomAccessGenotypeData, VariantFilterableGenotypeData {
	
	private final RandomAccessGenotypeData originalGenotypeData;
	private final VariantFilter variantFilter;

	public VariantFilterableGenotypeDataDecorator(RandomAccessGenotypeData originalGenotypeData, VariantFilter variantFilter) {
		
		if(variantFilter == null){
			throw new IllegalArgumentException("Variant filter can not be null");
		}
		
		this.originalGenotypeData = originalGenotypeData;
		this.variantFilter = variantFilter;
	}
	
	@Override
	public List<String> getSeqNames() {
		return originalGenotypeData.getSeqNames();
	}

	@Override
	public Iterable<Sequence> getSequences() {
		return originalGenotypeData.getSequences();
	}

	@Override
	public Sequence getSequenceByName(String name) {
		return originalGenotypeData.getSequenceByName(name);
	}

	@Override
	public Iterable<GeneticVariant> getVariantsByPos(String seqName, int startPos) {
		return createFilterIterable(originalGenotypeData.getVariantsByPos(seqName, startPos));
	}

	@Override
	public GeneticVariant getSnpVariantByPos(String seqName, int startPos) {
		GeneticVariant variant = originalGenotypeData.getSnpVariantByPos(seqName, startPos);
		
		if(variant == null){
			return null;
		}
		
		if(variantFilter.doesVariantPassFilter(variant)){
			return variant;
		} else {
			return null;
		}
	}

	@Override
	public Iterable<GeneticVariant> getSequenceGeneticVariants(String seqName) {
		return createFilterIterable(originalGenotypeData.getSequenceGeneticVariants(seqName));
	}

	@Override
	public Iterable<GeneticVariant> getVariantsByRange(String seqName, int rangeStart, int rangeEnd) {
		return createFilterIterable(originalGenotypeData.getVariantsByRange(seqName, rangeStart, rangeEnd));
	}

	@Override
	public List<Annotation> getVariantAnnotations() {
		return originalGenotypeData.getVariantAnnotations();
	}

	@Override
	public Annotation getVariantAnnotation(String annotationId) {
		return originalGenotypeData.getVariantAnnotation(annotationId);
	}

	@Override
	public List<SampleAnnotation> getSampleAnnotations() {
		return originalGenotypeData.getSampleAnnotations();
	}

	@Override
	public Annotation getSampleAnnotation(String annotationId) {
		return originalGenotypeData.getSampleAnnotation(annotationId);
	}

	@Override
	public List<Sample> getSamples() {
		return originalGenotypeData.getSamples();
	}

	@Override
	public Iterator<GeneticVariant> iterator() {
		return createFilterIterable(originalGenotypeData).iterator();
	}
	
	private Iterable<GeneticVariant> createFilterIterable(Iterable<GeneticVariant> originalIterable){
		return new VariantFilterIterable(new VariantFilterIterator(originalIterable.iterator(), variantFilter));
	}

	@Override
	public void setFilter(VariantFilter filter) {
		
		if(variantFilter == null){
			throw new IllegalArgumentException("Variant filter can not be null");
		}
	
	}

	@Override
	public VariantFilter getFilter() {
		return variantFilter;
	}
	
	@Override
	public void close() throws IOException {
		originalGenotypeData.close();
	}
	
}
