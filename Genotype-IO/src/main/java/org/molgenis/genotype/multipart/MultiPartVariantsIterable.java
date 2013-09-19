package org.molgenis.genotype.multipart;

import java.util.Iterator;
import java.util.NoSuchElementException;

import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;

public class MultiPartVariantsIterable implements Iterable<GeneticVariant>
{

	private final Iterable<RandomAccessGenotypeData> datasets;

	public MultiPartVariantsIterable(Iterable<RandomAccessGenotypeData> collection)
	{
		super();
		this.datasets = collection;

	}

	@Override
	public Iterator<GeneticVariant> iterator()
	{
		return new MultiPartVariantIterator(datasets);
	}

	private static class MultiPartVariantIterator implements Iterator<GeneticVariant>
	{

		private Iterator<RandomAccessGenotypeData> datasetIterator;
		private Iterator<GeneticVariant> datasetVariantIterator;
		private boolean done = false;

		public MultiPartVariantIterator(Iterable<RandomAccessGenotypeData> datasets)
		{
			super();
			this.datasetIterator = datasets.iterator();
			if (datasetIterator.hasNext())
			{
				this.datasetVariantIterator = datasetIterator.next().iterator();
			}
			else
			{
				done = true;
			}

		}

		@Override
		public boolean hasNext()
		{
			if (done)
			{
				return false;
			}

			while (!datasetVariantIterator.hasNext())
			{
				if (datasetIterator.hasNext())
				{
					datasetVariantIterator = datasetIterator.next().iterator();
				}
				else
				{
					done = true;
					return false;
				}
			}
			return true;

		}

		@Override
		public GeneticVariant next()
		{
			if (hasNext())
			{
				return datasetVariantIterator.next();
			}
			else
			{
				throw new NoSuchElementException();
			}
		}

		@Override
		public void remove()
		{
			throw new UnsupportedOperationException();
		}

	}

}
