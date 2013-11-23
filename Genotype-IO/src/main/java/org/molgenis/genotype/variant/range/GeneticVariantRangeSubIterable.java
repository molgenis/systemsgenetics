/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.variant.range;

import java.util.Iterator;
import java.util.List;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 *
 * @author Patrick Deelen
 */
public class GeneticVariantRangeSubIterable implements Iterable<GeneticVariant>{
	
	private final int listStartIndex;
	private final String seqStopName;
	private final int stopPos;
	private final List<GeneticVariant> sortedVariants;

	public GeneticVariantRangeSubIterable(int listStartIndex, String seqStopName, int stopPos, List<GeneticVariant> sortedVariants) {
		this.listStartIndex = listStartIndex;
		this.seqStopName = seqStopName;
		this.stopPos = stopPos;
		this.sortedVariants = sortedVariants;
	}

	@Override
	public Iterator<GeneticVariant> iterator() {
		return new GeneticVariantRangeSubIterator(listStartIndex, seqStopName, stopPos, sortedVariants);
	}
	
}
