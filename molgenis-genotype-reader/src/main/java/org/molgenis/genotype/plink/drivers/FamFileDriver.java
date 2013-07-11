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

import org.molgenis.genotype.plink.PlinkFileParser;
import org.molgenis.genotype.plink.datatypes.FamEntry;
import org.molgenis.util.TextFileUtils;

/**
 * Driver to query FAM files. FAM files annotate the families of BED files. See:
 * http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml
 * 
 * Content of a FAM file: Family ID Individual ID Paternal ID Maternal ID Sex
 * (1=male; 2=female; other=unknown) Phenotype
 */
public class FamFileDriver implements PlinkFileParser
{
	private BufferedReader reader;
	private File file;
	private String separators;
	private long nrElements;

	/**
	 * Construct a FamFileDriver on this file
	 * 
	 * @param famFile
	 * @throws Exception
	 */
	public FamFileDriver(File famFile)
	{
		this(famFile, DEFAULT_READ_FIELD_SEPARATORS);
	}

	public FamFileDriver(File famFile, char separator)
	{
		this(famFile, String.valueOf(separator));
	}
	
	public FamFileDriver(File famFile, String separators)
	{
		if (famFile == null) throw new IllegalArgumentException("file is null");
		this.file = famFile;
		this.separators = separators;
		this.nrElements = -1l;
	}


	/**
	 * Get all FAM file entries
	 * 
	 * @return
	 * @throws Exception
	 */
	public List<FamEntry> getAllEntries() throws IOException
	{
		reset();

		List<FamEntry> entryList = new ArrayList<FamEntry>();
		String line;
		while ((line = reader.readLine()) != null)
			entryList.add(parseEntry(line));

		return entryList;
	}

	/**
	 * Get a specific set of FAM file entries
	 * 
	 * @param from
	 *            = inclusive
	 * @param to
	 *            = exclusive
	 * @return
	 * @throws Exception
	 */
	public List<FamEntry> getEntries(final long from, final long to) throws IOException
	{
		reset();

		List<FamEntry> entryList = new ArrayList<FamEntry>();
		String line;
		for (int i = 0; (line = reader.readLine()) != null && i < to; ++i)
			if (i >= from) entryList.add(parseEntry(line));

		return entryList;
	}

	private FamEntry parseEntry(String line) throws IOException
	{
		StringTokenizer strTokenizer = new StringTokenizer(line, this.separators);
		try
		{
			String family = strTokenizer.nextToken();
			String individual = strTokenizer.nextToken();
			String father = strTokenizer.nextToken();
			String mother = strTokenizer.nextToken();
			byte sex = Byte.parseByte(strTokenizer.nextToken());
			double phenotype = Double.parseDouble(strTokenizer.nextToken());
			return new FamEntry(family, individual, father, mother, sex, phenotype);
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
