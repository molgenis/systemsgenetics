/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.sampleFilter;

import java.util.Iterator;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 *
 * @author Patrick Deelen
 */
public class ConvertToSampleFilterIterator implements Iterator<GeneticVariant> {

	private final Iterator<GeneticVariant> originalIterator;
	private final SampleFilterableGenotypeData genotypeData;

	public ConvertToSampleFilterIterator(Iterator<GeneticVariant> originalIterator, SampleFilterableGenotypeData genotypeData) {
		this.originalIterator = originalIterator;
		this.genotypeData = genotypeData;
	}

	@Override
	public boolean hasNext() {
		return originalIterator.hasNext();
	}

	@Override
	public GeneticVariant next() {
		return new SampleFilteredReadOnlyGeneticVariant(originalIterator.next(), genotypeData);
	}

	@Override
	public void remove() {
		throw new UnsupportedOperationException("Not supported. Use modifiable genotype data to remove individual variants or use gc genotype data");
	}
}
