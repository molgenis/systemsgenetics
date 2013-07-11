package org.molgenis.genotype.plink.writers;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;

import org.molgenis.genotype.plink.PlinkFileParser;
import org.molgenis.genotype.plink.datatypes.BimEntry;

/**
 * Write BIM file entries to a selected location.
 */
public class BimFileWriter implements PlinkFileParser
{
	private BufferedWriter writer;
	private char separator;

	public BimFileWriter(File bimFile) throws IOException
	{
		this(bimFile, DEFAULT_FIELD_SEPARATOR);
	}

	public BimFileWriter(File bimFile, char separator) throws IOException
	{
		if (bimFile == null) throw new IllegalArgumentException("file is null");
		this.writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(bimFile), FILE_ENCODING));
		this.separator = separator;
	}

	/**
	 * Write a single entry.
	 * 
	 * @throws IOException
	 */
	public void write(BimEntry bim) throws IOException
	{
		writer.write(bim.getChromosome());
		writer.write(this.separator);
		writer.write(bim.getSNP());
		writer.write(this.separator);
		writer.write(Double.toString(bim.getcM()));
		writer.write(this.separator);
		writer.write(Long.toString(bim.getBpPos()));
		writer.write(this.separator);
		writer.write(bim.getBiallele().get(0).toString());
		writer.write(this.separator);
		writer.write(bim.getBiallele().get(1).toString());
		writer.write(LINE_SEPARATOR);
	}

	/**
	 * Write multiple entries in order.
	 * 
	 * @throws IOException
	 */
	public void write(Iterable<BimEntry> bims) throws IOException
	{
		for (BimEntry bim : bims)
			write(bim);
	}

	/**
	 * Close the underlying writer.
	 */
	@Override
	public void close() throws IOException
	{
		if (writer != null) writer.close();
	}
}
