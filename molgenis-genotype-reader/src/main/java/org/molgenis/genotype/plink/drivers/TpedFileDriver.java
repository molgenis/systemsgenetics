package org.molgenis.genotype.plink.drivers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.StringTokenizer;

import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.plink.PlinkFileParser;
import org.molgenis.genotype.plink.datatypes.TpedEntry;
import org.molgenis.util.TextFileUtils;

/**
 * Driver to query TPED files.
 */
public class TpedFileDriver implements PlinkFileParser
{
	private BufferedReader reader;
	private File file;
	private String separators;
	private long nrElements;

	/**
	 * Construct a TpedFileDriver on this file
	 * 
	 * @param tpedFile
	 * @throws Exception
	 */
	public TpedFileDriver(File tpedFile)
	{
		this(tpedFile, DEFAULT_READ_FIELD_SEPARATORS);
	}

	public TpedFileDriver(File tpedFile, char separator)
	{
		this(tpedFile, String.valueOf(separator));
	}

	public TpedFileDriver(File tpedFile, String separators)
	{
		if (tpedFile == null) throw new IllegalArgumentException("file is null");
		this.file = tpedFile;
		this.separators = separators;
		this.nrElements = -1l;
	}

	/**
	 * Get a specific set of TPED file entries
	 * 
	 * @param from
	 *            = inclusive
	 * @param to
	 *            = exclusive
	 * @return
	 * @throws Exception
	 */
	public List<TpedEntry> getEntries(final long from, final long to) throws IOException
	{
		reset();

		List<TpedEntry> entryList = new ArrayList<TpedEntry>();
		String line;
		for (int i = 0; (line = reader.readLine()) != null && i < to; ++i)
			if (i >= from) entryList.add(parseEntry(line));

		return entryList;
	}

	/**
	 * Get all TPED file entries
	 * 
	 * @return
	 * @throws Exception
	 */
	public List<TpedEntry> getAllEntries() throws IOException
	{
		reset();

		List<TpedEntry> entryList = new ArrayList<TpedEntry>();
		String line;
		while ((line = reader.readLine()) != null)
			entryList.add(parseEntry(line));

		return entryList;
	}

	private TpedEntry parseEntry(String line) throws IOException
	{
		StringTokenizer strTokenizer = new StringTokenizer(line, separators);
		try
		{
			String chromosome = strTokenizer.nextToken();
			String snp = strTokenizer.nextToken();
			double cM = Double.parseDouble(strTokenizer.nextToken());
			long bpPos = Long.parseLong(strTokenizer.nextToken());
			List<Alleles> bialleles = new ArrayList<Alleles>();
			while (strTokenizer.hasMoreTokens())
			{
				char allele1 = strTokenizer.nextToken().charAt(0);
				char allele2 = strTokenizer.nextToken().charAt(0);
				bialleles.add(Alleles.createBasedOnChars(allele1, allele2));
			}
			return new TpedEntry(chromosome, snp, cM, bpPos, bialleles);
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
}
