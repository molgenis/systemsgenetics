package org.molgenis.genotype.multipart;

import java.util.Iterator;
import java.util.NoSuchElementException;

import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.Sequence;

public class MultiPartSequencesIterable implements Iterable<Sequence>
{

	private final Iterable<RandomAccessGenotypeData> datasets;

	public MultiPartSequencesIterable(Iterable<RandomAccessGenotypeData> collection)
	{
		super();
		this.datasets = collection;

	}

	@Override
	public Iterator<Sequence> iterator()
	{
		return new MultiPartSequencesIterator(datasets);
	}

	private class MultiPartSequencesIterator implements Iterator<Sequence>
	{

		private Iterator<RandomAccessGenotypeData> datasetIterator;
		private Iterator<Sequence> datasetSequenceIterator;
		private boolean done = false;

		public MultiPartSequencesIterator(Iterable<RandomAccessGenotypeData> datasets)
		{
			super();
			this.datasetIterator = datasets.iterator();
			if (datasetIterator.hasNext())
			{
				this.datasetSequenceIterator = datasetIterator.next().getSequences().iterator();
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

			while (!datasetSequenceIterator.hasNext())
			{
				if (datasetIterator.hasNext())
				{
					datasetSequenceIterator = datasetIterator.next().getSequences().iterator();
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
		public Sequence next()
		{
			if (hasNext())
			{
				return datasetSequenceIterator.next();
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
