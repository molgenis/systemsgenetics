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
import org.molgenis.genotype.plink.datatypes.BimEntry;
import org.molgenis.util.TextFileUtils;

/**
 * Driver to query BIM files. BIM files annotate the genotypes of BED files.
 * They are basically MAP files, with added biallelic data. See:
 * http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml
 * 
 * Content of a BIM file: chromosome, SNP, cM, base-position, allele 1, allele 2
 */
public class BimFileDriver implements PlinkFileParser
{
	private BufferedReader reader;
	private File file;
	private String separators;
	private long nrElements;

	/**
	 * Construct a BimFileDriver on this file
	 * 
	 * @param bimFile
	 * @throws Exception
	 */
	public BimFileDriver(File bimFile)
	{
		this(bimFile, DEFAULT_READ_FIELD_SEPARATORS);
	}

	public BimFileDriver(File bimFile, char separator)
	{
		this(bimFile, String.valueOf(separator));
	}

	public BimFileDriver(File bimFile, String separators)
	{
		if (bimFile == null) throw new IllegalArgumentException("file is null");
		this.file = bimFile;
		this.separators = separators;
		this.nrElements = -1l;
	}

	/**
	 * Get a specific set of BIM file entries
	 * 
	 * @param from
	 *            = inclusive
	 * @param to
	 *            = exclusive
	 * @return
	 * @throws Exception
	 */
	public List<BimEntry> getEntries(final long from, final long to) throws IOException
	{
		reset();

		List<BimEntry> entryList = new ArrayList<BimEntry>();
		String line;
		for (int i = 0; (line = reader.readLine()) != null && i < to; ++i)
			if (i >= from) entryList.add(parseEntry(line));

		return entryList;
	}

	/**
	 * Get all BIM file entries
	 * 
	 * @return
	 * @throws Exception
	 */
	public List<BimEntry> getAllEntries() throws IOException
	{
		reset();

		List<BimEntry> entryList = new ArrayList<BimEntry>();
		String line;
		while ((line = reader.readLine()) != null)
			entryList.add(parseEntry(line));

		return entryList;
	}

	private BimEntry parseEntry(String line) throws IOException
	{
		StringTokenizer strTokenizer = new StringTokenizer(line, this.separators);
		try
		{
			String chromosome = strTokenizer.nextToken();
			String snp = strTokenizer.nextToken();
			double cM = Double.parseDouble(strTokenizer.nextToken());
			long bpPos = Long.parseLong(strTokenizer.nextToken());
			char allelle1 = strTokenizer.nextToken().charAt(0);
			char allelle2 = strTokenizer.nextToken().charAt(0);
			return new BimEntry(chromosome, snp, cM, bpPos, Alleles.createBasedOnChars(allelle1, allelle2));
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
