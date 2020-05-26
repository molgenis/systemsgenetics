package org.molgenis.genotype.tabix;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

import net.sf.samtools.util.BlockCompressedInputStream;

import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.RawLineQuery;
import org.molgenis.genotype.RawLineQueryResult;
import org.molgenis.genotype.tabix.TabixIndex.TabixIterator;

/**
 * Returns raw lines instead of GeneticVariants
 * 
 * @author erwin
 * 
 */
public class TabixRawLineQuery implements RawLineQuery
{
	private BlockCompressedInputStream inputStream;
	private final TabixIndex index;

	public TabixRawLineQuery(File bzipFile, TabixIndex index)
	{
		if (bzipFile == null) throw new IllegalArgumentException("BzipFile is null");
		if (index == null) throw new IllegalArgumentException("Index is null");

		try
		{
			inputStream = new BlockCompressedInputStream(bzipFile);
		}
		catch (IOException e)
		{
			throw new GenotypeDataException("IOExcption creating BlockCompressedInputStream for " + bzipFile.getName(),
					e);
		}

		this.index = index;
	}

	@Override
	public RawLineQueryResult executeQuery(String sequence, int startPos)
	{
		if (startPos < 0) throw new IllegalArgumentException("StartPos must be bigger then 0");

		try
		{
			TabixIterator tabixIterator = index.queryTabixIndex(sequence, startPos - 1, startPos, inputStream);
			return new TabixRawLineQueryResult(inputStream, new RawLineIterator(tabixIterator));
		}
		catch (IOException e)
		{
			throw new GenotypeDataException(e);
		}
	}

	private static class RawLineIterator implements Iterator<String>
	{
		private TabixIterator tabixIterator;
		private String line;

		private RawLineIterator(TabixIterator tabixIterator) throws IOException
		{
			super();
			this.tabixIterator = tabixIterator;
			line = tabixIterator == null ? null : tabixIterator.next();
		}

		@Override
		public boolean hasNext()
		{
			return line != null;
		}

		@Override
		public String next()
		{
			String result = line;

			try
			{
				line = tabixIterator.next();
			}
			catch (IOException e)
			{
				throw new RuntimeException("Exception calling next on TabixIndex.TabixIterator", e);
			}

			return result;
		}

		@Override
		public void remove()
		{
			throw new UnsupportedOperationException();
		}

	}
}
