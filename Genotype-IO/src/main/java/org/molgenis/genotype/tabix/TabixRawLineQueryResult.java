package org.molgenis.genotype.tabix;

import java.io.IOException;
import java.io.InputStream;
import java.util.Iterator;

import org.molgenis.genotype.RawLineQueryResult;

public class TabixRawLineQueryResult implements RawLineQueryResult
{
	private InputStream inputStream;
	private Iterator<String> iterator;

	public TabixRawLineQueryResult(InputStream inputStream, Iterator<String> iterator)
	{
		super();
		this.inputStream = inputStream;
		this.iterator = iterator;
	}

	@Override
	public void close() throws IOException
	{
		inputStream.close();

	}

	@Override
	public Iterator<String> iterator()
	{
		return iterator;
	}

}
