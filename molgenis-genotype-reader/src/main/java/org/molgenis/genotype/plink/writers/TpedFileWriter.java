package org.molgenis.genotype.plink.writers;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;

import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.plink.PlinkFileParser;
import org.molgenis.genotype.plink.datatypes.TpedEntry;

/**
 * Write MAP file entries to a selected location.
 */
public class TpedFileWriter implements PlinkFileParser
{
	private BufferedWriter writer;
	private char separator;

	public TpedFileWriter(File tpedFile) throws IOException
	{
		this(tpedFile, DEFAULT_FIELD_SEPARATOR);
	}

	public TpedFileWriter(File tpedFile, char separator) throws IOException
	{
		if (tpedFile == null) throw new IllegalArgumentException("tped file is null");
		this.writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(tpedFile), FILE_ENCODING));
		this.separator = separator;
	}

	/**
	 * Write a single entry.
	 * 
	 * @throws IOException
	 */
	public void write(TpedEntry tped) throws IOException
	{
		writer.write(tped.getChromosome());
		writer.write(separator);
		writer.write(tped.getSNP());
		writer.write(separator);
		writer.write(Double.toString(tped.getcM()));
		writer.write(separator);
		writer.write(Long.toString(tped.getBpPos()));
		for (Alleles biallele : tped.getBialleles())
		{
			writer.write(separator);
			writer.write(biallele.get(0).toString());
			writer.write(separator);
			writer.write(biallele.get(1).toString());
		}
		writer.write(LINE_SEPARATOR);
	}

	/**
	 * Write multiple entries in order.
	 * 
	 * @throws IOException
	 */
	public void write(Iterable<TpedEntry> tpeds) throws IOException
	{
		for (TpedEntry tped : tpeds)
			write(tped);
	}

	/**
	 * Close the underlying writer.
	 * 
	 * @throws IOException
	 */
	@Override
	public void close() throws IOException
	{
		if (writer != null) writer.close();
	}
}
