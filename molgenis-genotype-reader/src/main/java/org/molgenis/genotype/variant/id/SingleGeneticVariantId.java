package org.molgenis.genotype.variant.id;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

public class SingleGeneticVariantId extends GeneticVariantId
{

	private final String variantId;

	protected SingleGeneticVariantId(String variantId)
	{
		super();
		this.variantId = variantId;
	}

	@Override
	public String getPrimairyId()
	{
		return variantId;
	}

	@Override
	public List<String> getVariantIds()
	{
		ArrayList<String> variantIds = new ArrayList<String>(1);
		variantIds.add(variantId);
		return Collections.unmodifiableList(variantIds);
	}

	@Override
	public String getConcatenatedId()
	{
		return variantId;
	}

	@Override
	public String getConcatenatedId(String separator)
	{
		return variantId;
	}

	@Override
	public boolean isIdInVariantIds(String queryId)
	{
		return variantId.equals(queryId);
	}

	@Override
	public boolean onlyPrimairyId()
	{
		return true;
	}

	@Override
	public Iterator<String> iterator()
	{
		return getVariantIds().iterator();
	}

	@Override
	public List<String> getAlternativeIds()
	{
		return Collections.emptyList();
	}

	@Override
	public boolean containsId()
	{
		return true;
	}

	@Override
	public int hashCode()
	{
		final int prime = 31;
		int result = 1;
		result = prime * result + ((variantId == null) ? 0 : variantId.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj)
	{
		if (this == obj) return true;
		if (obj == null) return false;
		if (getClass() != obj.getClass()) return false;
		SingleGeneticVariantId other = (SingleGeneticVariantId) obj;
		if (variantId == null)
		{
			if (other.variantId != null) return false;
		}
		else if (!variantId.equals(other.variantId)) return false;
		return true;
	}

}
