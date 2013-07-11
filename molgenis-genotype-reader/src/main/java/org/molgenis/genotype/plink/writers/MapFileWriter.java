package org.molgenis.genotype.plink.writers;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;

import org.molgenis.genotype.plink.PlinkFileParser;
import org.molgenis.genotype.plink.datatypes.MapEntry;

/**
 * Write MAP file entries to a selected location.
 */
public class MapFileWriter implements PlinkFileParser
{
	private BufferedWriter writer;
	private char separator;

	public MapFileWriter(File mapFile) throws IOException
	{
		this(mapFile, DEFAULT_FIELD_SEPARATOR);
	}

	public MapFileWriter(File mapFile, char separator) throws IOException
	{
		if (mapFile == null) throw new IllegalArgumentException("file is null");
		this.writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(mapFile), FILE_ENCODING));
		this.separator = separator;
	}

	/**
	 * Write a single entry.
	 * 
	 * @throws IOException
	 */
	public void write(MapEntry map) throws IOException
	{
		writer.write(map.getChromosome());
		writer.write(separator);
		writer.write(map.getSNP());
		writer.write(separator);
		if (map.getcM() == 0d)
		{
			writer.write("0");
		}
		else
		{
			writer.write(Double.toString(map.getcM()));
		}
		writer.write(separator);
		writer.write(Long.toString(map.getBpPos()));
		writer.write(LINE_SEPARATOR);
	}

	/**
	 * Write multiple entries in order.
	 * 
	 * @throws IOException
	 */
	public void write(Iterable<MapEntry> maps) throws IOException
	{
		for (MapEntry map : maps)
			write(map);
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
