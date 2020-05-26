package org.molgenis.genotype.plink.readers;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.util.Iterator;

import org.molgenis.genotype.plink.PlinkFileParser;
import org.molgenis.genotype.plink.datatypes.MapEntry;
import org.molgenis.genotype.plink.drivers.MapFileDriver;

public class MapFileReader implements PlinkFileParser, Iterable<MapEntry>
{
	private static final Charset CHARSET_UTF8 = Charset.forName("UTF-8");
	private final BufferedReader reader;
	private final String separators;

	public MapFileReader(InputStream in)
	{
		this(in, DEFAULT_READ_FIELD_SEPARATORS);
	}
	
	public MapFileReader(InputStream in, char separator)
	{
		this(in, String.valueOf(separator));
	}
	
	public MapFileReader(InputStream in, String separators)
	{
		reader = new BufferedReader(new InputStreamReader(in, CHARSET_UTF8));
		this.separators = separators;
	}

	@Override
	public void close() throws IOException
	{
		reader.close();
	}

	@Override
	public Iterator<MapEntry> iterator()
	{
		return new MapFileIterator();
	}

	private class MapFileIterator implements Iterator<MapEntry>
	{
		private String line;

		public MapFileIterator()
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
		public MapEntry next()
		{
			MapEntry entry;
			try
			{
				entry = MapFileDriver.parseEntry(line, separators);
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
