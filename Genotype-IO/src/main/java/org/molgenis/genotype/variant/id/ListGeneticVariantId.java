package org.molgenis.genotype.variant.id;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

public class ListGeneticVariantId extends GeneticVariantId
{

	private final List<String> variantIds;

	protected ListGeneticVariantId(List<String> variantIds)
	{
		super();

		if (variantIds.size() < 2)
		{
			throw new IllegalArgumentException("Use SingleGeneticVariantId in case of single variant ID");
		}

		this.variantIds = variantIds;
	}

	protected ListGeneticVariantId(String primairyId, List<String> alternativeIds)
	{

		if (alternativeIds == null)
		{
			throw new IllegalArgumentException("Alternative IDs list is null");
		}
		if (alternativeIds.isEmpty())
		{
			throw new IllegalArgumentException("Alternative IDs list is empty");
		}

		ArrayList<String> variantIdsBuilder = new ArrayList<String>(alternativeIds.size() + 1);
		variantIdsBuilder.add(primairyId);
		variantIdsBuilder.addAll(alternativeIds);

		this.variantIds = variantIdsBuilder;
	}

	@Override
	public String getPrimairyId()
	{
		return variantIds.get(0);
	}

	@Override
	public List<String> getVariantIds()
	{
		return Collections.unmodifiableList(variantIds);
	}

	@Override
	public String getConcatenatedId()
	{
		return getConcatenatedId(DEFAULT_ID_SEPARATOR);
	}

	@Override
	public String getConcatenatedId(String separator)
	{
		StringBuilder concatenatedIds = new StringBuilder();

		boolean notFirst = false;
		for (String variantId : variantIds)
		{
			if (notFirst)
			{
				concatenatedIds.append(separator);
			}
			concatenatedIds.append(variantId);
		}

		return concatenatedIds.toString();

	}

	@Override
	public boolean isIdInVariantIds(String queryId)
	{
		for (String variantId : this)
		{
			if (variantId.equals(queryId))
			{
				return true;
			}
		}
		return false;
	}

	@Override
	public boolean onlyPrimairyId()
	{
		return false;
	}

	@Override
	public Iterator<String> iterator()
	{
		return getVariantIds().iterator();
	}

	@Override
	public List<String> getAlternativeIds()
	{
		return Collections.unmodifiableList(variantIds.subList(1, variantIds.size()));
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
		result = prime * result + ((variantIds == null) ? 0 : variantIds.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj)
	{
		if (this == obj) return true;
		if (obj == null) return false;
		if (getClass() != obj.getClass()) return false;
		ListGeneticVariantId other = (ListGeneticVariantId) obj;
		if (variantIds == null)
		{
			if (other.variantIds != null) return false;
		}
		else if (!variantIds.equals(other.variantIds)) return false;
		return true;
	}

}
