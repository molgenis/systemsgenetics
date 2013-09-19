package org.molgenis.genotype.tabix;

import java.io.IOException;
import java.io.InputStream;
import java.util.Iterator;

import org.molgenis.genotype.VariantQueryResult;
import org.molgenis.genotype.variant.GeneticVariant;

public class TabixQueryResult implements VariantQueryResult
{
	private final InputStream inputStream;
	private final Iterator<GeneticVariant> iterator;

	public TabixQueryResult(InputStream inputStream, Iterator<GeneticVariant> iterator)
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
	public Iterator<GeneticVariant> iterator()
	{
		return iterator;
	}

}
