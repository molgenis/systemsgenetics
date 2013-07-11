package org.molgenis.genotype.plink.writers;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;

import org.molgenis.genotype.plink.PlinkFileParser;
import org.molgenis.genotype.plink.datatypes.FamEntry;

/**
 * Write MAP file entries to a selected location.
 */
public class FamFileWriter implements PlinkFileParser
{
	private BufferedWriter writer;
	private char separator;

	public FamFileWriter(File famFile) throws IOException
	{
		this(famFile, DEFAULT_FIELD_SEPARATOR);
	}

	public FamFileWriter(File famFile, char separator) throws IOException
	{
		if (famFile == null) throw new IllegalArgumentException("file is null");
		this.writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(famFile), FILE_ENCODING));
		this.separator = separator;
	}

	/**
	 * Write a single entry.
	 * 
	 * @throws IOException
	 */
	public void write(FamEntry fam) throws IOException
	{
		writer.write(fam.getFamily());
		writer.write(separator);
		writer.write(fam.getIndividual());
		writer.write(separator);
		writer.write(fam.getFather());
		writer.write(separator);
		writer.write(fam.getMother());
		writer.write(separator);
		writer.write(Byte.toString(fam.getSex()));
		writer.write(separator);
		writer.write(Double.toString(fam.getPhenotype()));
		writer.write(LINE_SEPARATOR);
	}

	/**
	 * Write multiple entries in order.
	 * 
	 * @throws IOException
	 */
	public void write(Iterable<FamEntry> fams) throws IOException
	{
		for (FamEntry fam : fams)
			write(fam);
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
