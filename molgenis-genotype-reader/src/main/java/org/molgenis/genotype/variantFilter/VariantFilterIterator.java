/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.variantFilter;

import java.util.Iterator;
import java.util.NoSuchElementException;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 *
 * @author Patrick Deelen
 */
public class VariantFilterIterator implements Iterator<GeneticVariant> {
	
	private final Iterator<GeneticVariant> originalIterator;
	private final VariantFilter qcChecker;
			
	private GeneticVariant next;

	public VariantFilterIterator(Iterator<GeneticVariant> originalIterator, VariantFilter qcChecker) {
		this.originalIterator = originalIterator;
		this.qcChecker = qcChecker;
	}
	
	@Override
	public boolean hasNext()
	{
		return next != null;
	}
	
	@Override
	public GeneticVariant next() {
		
		if (next == null)
		{
			throw new NoSuchElementException();
		}

		GeneticVariant currentNext = next;

		// prepare next next
		goToNext();

		return currentNext;
	}
	
	private void goToNext()
	{
		while (originalIterator.hasNext())
		{
			GeneticVariant originalNext = originalIterator.next();

			if (!qcChecker.doesVariantPassFilter(originalNext))
			{
				// skip variants on exclude list
				continue;
			}
			
			next = originalNext;
			return;
		}
		// We do a return if we find a non excluded next. So if we get here it
		// is the end of the original iterator. Setting next to null so hasNext
		// knows it is the end.
		next = null;
	}



	@Override
	public void remove() {
		throw new UnsupportedOperationException("Not supported yet.");
	}
	
}
