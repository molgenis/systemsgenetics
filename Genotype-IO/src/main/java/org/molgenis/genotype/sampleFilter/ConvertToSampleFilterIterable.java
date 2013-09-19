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
public class ConvertToSampleFilterIterable implements Iterable<GeneticVariant> {

	private final Iterable<GeneticVariant> originalIterable;
	private final SampleFilterableGenotypeData genotypeData;

	public ConvertToSampleFilterIterable(Iterable<GeneticVariant> originalIterable, SampleFilterableGenotypeData genotypeData) {
		this.originalIterable = originalIterable;
		this.genotypeData = genotypeData;
	}

	@Override
	public Iterator<GeneticVariant> iterator() {
		return new ConvertToSampleFilterIterator(originalIterable.iterator(), genotypeData);
	}
}
