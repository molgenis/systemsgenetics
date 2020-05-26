package org.molgenis.genotype.variant.id;

import java.util.Collections;
import java.util.Iterator;
import java.util.List;

public class BlankGeneticVariantId extends GeneticVariantId
{

	protected BlankGeneticVariantId()
	{
	}

	@Override
	public Iterator<String> iterator()
	{
		return this.getVariantIds().iterator();
	}

	@Override
	public String getPrimairyId()
	{
		return null;
	}

	@Override
	public List<String> getVariantIds()
	{
		return Collections.emptyList();
	}

	@Override
	public List<String> getAlternativeIds()
	{
		return Collections.emptyList();
	}

	@Override
	public String getConcatenatedId()
	{
		return "";
	}

	@Override
	public String getConcatenatedId(String separator)
	{
		return "";
	}

	@Override
	public boolean isIdInVariantIds(String queryId)
	{
		return false;
	}

	@Override
	public boolean onlyPrimairyId()
	{
		return false;
	}

	@Override
	public boolean containsId()
	{
		return false;
	}

	@Override
	public int hashCode()
	{
		return 1;
	}

	@Override
	public boolean equals(Object obj)
	{
		if (this == obj) return true;
		if (obj == null) return false;
		if (getClass() != obj.getClass()) return false;

		return true;
	}

}
