package org.molgenis.genotype.modifiable;

import java.util.HashSet;
import java.util.Iterator;
import java.util.NoSuchElementException;

import org.molgenis.genotype.variant.GeneticVariant;

public class ModifiableGeneticVariantIterator<E extends GeneticVariant> implements Iterator<E>
{

	private final Iterator<GeneticVariant> originalIterator;
	private final ModifiableGenotypeData modifiableGenotypeData;
	private final HashSet<ModifiableGeneticVariant> excludeList;
	private ModifiableGeneticVariant next;

	public ModifiableGeneticVariantIterator(Iterator<GeneticVariant> originalIterator,
			ModifiableGenotypeData modifiableGenotypeData, HashSet<ModifiableGeneticVariant> excludeList)
	{
		super();
		this.originalIterator = originalIterator;
		this.modifiableGenotypeData = modifiableGenotypeData;
		this.excludeList = excludeList;

		goToNext();
	}

	@Override
	public boolean hasNext()
	{
		return next != null;
	}

	@SuppressWarnings("unchecked")
	@Override
	public E next()
	{
		if (next == null)
		{
			throw new NoSuchElementException();
		}

		ModifiableGeneticVariant currentNext = next;

		// prepare next next
		goToNext();

		return (E) currentNext;

	}

	@Override
	public void remove()
	{
		throw new UnsupportedOperationException("Removal not supported via iterator.");
	}

	private void goToNext()
	{
		while (originalIterator.hasNext())
		{
			ModifiableGeneticVariant originalNext = new ModifiableGeneticVariant(originalIterator.next(),
					modifiableGenotypeData);

			if (excludeList.contains(originalNext))
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

	/**
	 * Wrap genetic variant iterator to return modifiable genetic variants
	 * 
	 * @param originalIterator
	 *            the original iterator
	 * @param modifiableGenotypeData
	 *            the modifiable genotype data that stores the changes.
	 * @param excludeList
	 *            variant to skip over when iterating
	 * @return
	 */
	public static Iterable<ModifiableGeneticVariant> createModifiableGeneticVariantIterable(
			Iterator<GeneticVariant> originalIterator, ModifiableGenotypeData modifiableGenotypeData,
			HashSet<ModifiableGeneticVariant> excludeList)
	{
		return new ModifiableGeneticVariantIterable<ModifiableGeneticVariant>(
				new ModifiableGeneticVariantIterator<ModifiableGeneticVariant>(originalIterator,
						modifiableGenotypeData, excludeList));
	}

	/**
	 * Wrap genetic variant iterable to return genetic variants that are
	 * modifiable
	 * 
	 * @param originalIterator
	 * @param modifiableGenotypeData
	 * @param excludeList
	 *            variant to skip over when iterating
	 * @return
	 */
	public static Iterable<GeneticVariant> createGeneticVariantIterableBackByModifiable(
			Iterator<GeneticVariant> originalIterator, ModifiableGenotypeData modifiableGenotypeData,
			HashSet<ModifiableGeneticVariant> excludeList)
	{
		return new ModifiableGeneticVariantIterable<GeneticVariant>(
				new ModifiableGeneticVariantIterator<GeneticVariant>(originalIterator, modifiableGenotypeData,
						excludeList));
	}

	protected static class ModifiableGeneticVariantIterable<E extends GeneticVariant> implements Iterable<E>
	{

		private final Iterator<E> modifiableGeneticVariantIterator;

		public ModifiableGeneticVariantIterable(Iterator<E> modifiableGeneticVariantIterator)
		{
			super();
			this.modifiableGeneticVariantIterator = modifiableGeneticVariantIterator;
		}

		@Override
		public Iterator<E> iterator()
		{
			return modifiableGeneticVariantIterator;
		}

	}
}
