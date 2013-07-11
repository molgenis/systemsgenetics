package org.molgenis.genotype.plink.writers;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;

import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.plink.PlinkFileParser;
import org.molgenis.genotype.plink.datatypes.PedEntry;

/**
 * Write MAP file entries to a selected location.
 */
public class PedFileWriter implements PlinkFileParser
{
	private final BufferedWriter writer;
	private final char separator;

	public PedFileWriter(File pedFile) throws IOException
	{
		this(pedFile, DEFAULT_FIELD_SEPARATOR);
	}

	public PedFileWriter(File pedFile, char separator) throws IOException
	{
		if (pedFile == null) throw new IllegalArgumentException("ped file is null");
		this.writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(pedFile), FILE_ENCODING));
		this.separator = separator;
	}

	/**
	 * Write a single entry.
	 * 
	 * @throws IOException
	 */
	public void write(PedEntry ped) throws IOException
	{
		writer.write(ped.getFamily());
		writer.write(separator);
		writer.write(ped.getIndividual());
		writer.write(separator);
		writer.write(ped.getFather());
		writer.write(separator);
		writer.write(ped.getMother());
		writer.write(separator);
		writer.write(Byte.toString(ped.getSex()));
		writer.write(separator);
		writer.write(Double.toString(ped.getPhenotype()));

		for (Alleles biallele : ped)
		{
			writer.write(separator);
			writer.write(biallele.get(0).toString());
			writer.write(separator);
			writer.write(biallele.get(1).toString());
		}

		writer.write(LINE_SEPARATOR);
		writer.flush();
	}

	/**
	 * Write multiple entries in order.
	 * 
	 * @throws IOException
	 */
	public void write(Iterable<PedEntry> peds) throws IOException
	{
		for (PedEntry ped : peds)
			write(ped);
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
