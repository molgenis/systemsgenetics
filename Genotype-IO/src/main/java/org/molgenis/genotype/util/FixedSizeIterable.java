package org.molgenis.genotype.util;

/**
 * Iterable that iterates over a fixed number of elements
 */
public interface FixedSizeIterable<E> extends Iterable<E>
{
	/**
	 * @return the number of items for this iterable
	 */
	public int size();
}
