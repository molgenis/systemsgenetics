package org.molgenis.genotype.plink.drivers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.StringTokenizer;

import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.plink.PlinkFileParser;
import org.molgenis.genotype.plink.datatypes.PedEntry;
import org.molgenis.util.TextFileUtils;

/**
 * Driver to query PED files. A PED file contains family- and genotyping data
 * for an individual, plus a single phenotype. Basically it is a FAM file with
 * added genotyping (typically SNP) data. However, the example file is a bit
 * peculiar: it has 'null' columns because of additional spacing between some
 * data values. This makes parsing hard. Question: can all Plink files have
 * this? or just PED? See:
 * http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped
 * 
 * Answer by Patrick: This (should) never really happens. Only if you copy paste
 * from the website where they use spaces to make it look nice.
 */
public class PedFileDriver implements PlinkFileParser, Iterable<PedEntry>
{
	private BufferedReader reader;
	private final File file;
	private final String separators;
	private long nrElements;

	/**
	 * Construct a PedFileDriver on this file
	 * 
	 * @param bimFile
	 * @throws Exception
	 */
	public PedFileDriver(File pedFile)
	{
		this(pedFile, DEFAULT_READ_FIELD_SEPARATORS);
	}

	public PedFileDriver(File pedFile, char separator)
	{
		this(pedFile, String.valueOf(separator));
	}

	public PedFileDriver(File pedFile, String separators)
	{
		if (pedFile == null) throw new IllegalArgumentException("file is null");
		this.file = pedFile;
		this.separators = separators;
		this.nrElements = -1l;
	}

	/**
	 * Get all PED file entries
	 * 
	 * @return
	 * @throws Exception
	 */
	public List<PedEntry> getAllEntries() throws IOException
	{
		reset();

		List<PedEntry> entryList = new ArrayList<PedEntry>();
		String line;
		while ((line = reader.readLine()) != null)
			entryList.add(parseEntry(line));

		return entryList;
	}

	@Override
	public Iterator<PedEntry> iterator()
	{
		try
		{
			reset();
			return new PedFileIterator();
		}
		catch (IOException e)
		{
			throw new RuntimeException(e);
		}
	}

	/**
	 * Get a specific set of PED file entries
	 * 
	 * @param from
	 *            = inclusive
	 * @param to
	 *            = exclusive
	 * @return
	 * @throws Exception
	 */
	public List<PedEntry> getEntries(final long from, final long to) throws IOException
	{
		reset();

		List<PedEntry> entryList = new ArrayList<PedEntry>();
		String line;
		for (int i = 0; (line = reader.readLine()) != null && i < to; ++i)
			if (i >= from) entryList.add(parseEntry(line));

		return entryList;
	}

	private PedEntry parseEntry(String line) throws IOException
	{
		final StringTokenizer strTokenizer = new StringTokenizer(line, separators);
		try
		{
			//Todo this to make sure not to save whole line in background after a split or tokanizer
			String family = new String(strTokenizer.nextToken());
			String individual = new String(strTokenizer.nextToken());
			String father = new String(strTokenizer.nextToken());
			String mother = new String(strTokenizer.nextToken());
			byte sex = Byte.parseByte(strTokenizer.nextToken());
			double phenotype = Double.parseDouble(strTokenizer.nextToken());

			return new PedEntry(family, individual, father, mother, sex, phenotype, new Iterator<Alleles>()
			{

				@Override
				public boolean hasNext()
				{
					return strTokenizer.hasMoreTokens();
				}

				@Override
				public Alleles next()
				{
					char allele1 = strTokenizer.nextToken().charAt(0);
					char allele2 = strTokenizer.nextToken().charAt(0);
					return Alleles.createBasedOnChars(allele1, allele2);
				}

				@Override
				public void remove()
				{
					throw new UnsupportedOperationException();
				}

			});
		}
		catch (NoSuchElementException e)
		{
			throw new IOException("error in line: " + line, e);
		}
		catch (IndexOutOfBoundsException e)
		{
			throw new IOException("error in line: " + line, e);
		}
		catch (NumberFormatException e)
		{
			throw new IOException("error in line: " + line, e);
		}
	}

	public long getNrOfElements() throws IOException
	{
		if (nrElements == -1) nrElements = TextFileUtils.getNumberOfNonEmptyLines(file, FILE_ENCODING);
		return nrElements;
	}

	@Override
	public void close() throws IOException
	{
		if (this.reader != null) this.reader.close();
	}

	public void reset() throws IOException
	{
		if (this.reader != null) close();
		this.reader = new BufferedReader(new InputStreamReader(new FileInputStream(file), FILE_ENCODING));
	}

	private class PedFileIterator implements Iterator<PedEntry>
	{
		private String line;

		public PedFileIterator()
		{
			try
			{
				line = reader.readLine();
			}
			catch (IOException e)
			{
				throw new RuntimeException(e);
			}
		}

		@Override
		public boolean hasNext()
		{
			return line != null;
		}

		@Override
		public PedEntry next()
		{
			PedEntry entry;
			try
			{
				entry = parseEntry(line);
				line = reader.readLine();
			}
			catch (IOException e)
			{
				throw new RuntimeException(e);
			}

			return entry;
		}

		@Override
		public void remove()
		{
			throw new UnsupportedOperationException();
		}

	}
}
