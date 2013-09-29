/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.variantFilter;

import java.io.IOException;
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
 * View of orginal genotype data only including variants conforming to set
 * filter values
 *
 *
 * @author Patrick Deelen
 */
public class RandomAccessGenotypeDataVariantQc extends AbstractRandomAccessGenotypeData {

	private final RandomAccessGenotypeData originalGenotypeData;
	private final VariantQcChecker qcChecker;

	public RandomAccessGenotypeDataVariantQc(RandomAccessGenotypeData originalGenotypeData, VariantQcChecker qcChecker) {
		this.originalGenotypeData = originalGenotypeData;
		this.qcChecker = qcChecker;
	}
	
	public RandomAccessGenotypeDataVariantQc(RandomAccessGenotypeData originalGenotypeData){
		this.originalGenotypeData = originalGenotypeData;
		this.qcChecker = new VariantQcChecker(1, 0, 0);
	}
	
	public void setMafCutoff(float maf){
		this.qcChecker.setMafCutoff(maf);
	}

	public void setHweCutoff(double hwe) {
		this.qcChecker.setHweCutoff(hwe);
	}

	public void setCallRateCutoff(double callRate) {
		this.qcChecker.setCallRateCutoff(callRate);
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
		return createQcIterable(originalGenotypeData.getVariantsByPos(seqName, startPos));
	}

	@Override
	public GeneticVariant getSnpVariantByPos(String seqName, int startPos) {
		GeneticVariant variant = originalGenotypeData.getSnpVariantByPos(seqName, startPos);
		
		if(variant == null){
			return null;
		}
		
		if(qcChecker.doesVariantPassFilter(variant)){
			return variant;
		} else {
			return null;
		}
	}

	@Override
	public Iterable<GeneticVariant> getSequenceGeneticVariants(String seqName) {
		return createQcIterable(originalGenotypeData.getSequenceGeneticVariants(seqName));
	}

	@Override
	public Iterable<GeneticVariant> getVariantsByRange(String seqName, int rangeStart, int rangeEnd) {
		return createQcIterable(originalGenotypeData.getVariantsByRange(seqName, rangeStart, rangeEnd));
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
		return createQcIterable(originalGenotypeData).iterator();
	}
	
	private Iterable<GeneticVariant> createQcIterable(Iterable<GeneticVariant> originalIterable){
		return new VariantFilterIterable(new VariantFilterIterator(originalIterable.iterator(), qcChecker));
	}

	@Override
	public void close() throws IOException {
	}

	@Override
	public Map<String, Annotation> getVariantAnnotationsMap() {
		return originalGenotypeData.getVariantAnnotationsMap();
	}

	@Override
	public Map<String, SampleAnnotation> getSampleAnnotationsMap() {
		return originalGenotypeData.getSampleAnnotationsMap();
	}

	@Override
	public boolean isOnlyContaingSaveProbabilityGenotypes() {
		return originalGenotypeData.isOnlyContaingSaveProbabilityGenotypes();
	}
	
}
