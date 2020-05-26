/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.variant.range;

import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import org.molgenis.genotype.util.ChromosomeComparator;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 *
 * @author Patrick Deelen
 */
public class GeneticVariantRangeSubIterator implements Iterator<GeneticVariant>{

	private int currentIndex;
	private final String seqStopName;
	private final int stopPos;
	private final List<GeneticVariant> sortedVariants;

	public GeneticVariantRangeSubIterator(int listStartIndex, String seqStopName, int stopPos, List<GeneticVariant> sortedVariants) {
		this.currentIndex = listStartIndex - 1; // do -1 so next will include the first index of interest
		this.seqStopName = seqStopName;
		this.stopPos = stopPos;
		this.sortedVariants = sortedVariants;		
	}
	
	@Override
	public boolean hasNext() {
		
		if(currentIndex + 1 < sortedVariants.size()){
			GeneticVariant nextVariant = sortedVariants.get(currentIndex + 1);
			
			if (nextVariant.getSequenceName().equals(seqStopName) && nextVariant.getStartPos() <= stopPos ){
				return true;
			}  else if(ChromosomeComparator.chrASmallerChrB(nextVariant.getSequenceName(), seqStopName)){
				return true;
			} else {
				return false;
			}
			
		} else {
			return false;
		}
		
	}

	@Override
	public GeneticVariant next() {
		if(!hasNext()){
			throw new NoSuchElementException("No more variants in selected range");
		} else {
			++currentIndex;
			return sortedVariants.get(currentIndex);
		}
	}

	@Override
	public void remove() {
		throw new UnsupportedOperationException("Not supported ever."); 
	}

}
