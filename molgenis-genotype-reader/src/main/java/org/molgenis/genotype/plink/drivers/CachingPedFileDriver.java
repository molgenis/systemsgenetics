package org.molgenis.genotype.plink.drivers;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.molgenis.framework.db.QueryRule;
import org.molgenis.framework.db.QueryRule.Operator;
import org.molgenis.genotype.plink.datatypes.PedEntry;

public class CachingPedFileDriver extends PedFileDriver
{
	private List<PedEntry> entries;
	private List<PedEntry> filteredEntries;

	public CachingPedFileDriver(File pedFile)
	{
		super(pedFile, DEFAULT_FIELD_SEPARATOR);
	}

	public CachingPedFileDriver(File pedFile, char separator)
	{
		super(pedFile, separator);
	}

	public void setFilters(List<QueryRule> rules, List<String> snpNames) throws IOException
	{
		// lazy initialization
		if (entries == null) entries = super.getAllEntries();
		this.filteredEntries = new ArrayList<PedEntry>(entries);

		for (QueryRule rule : rules)
		{
			filteredEntries = filter(filteredEntries, rule, snpNames);
		}
	}

	@Override
	public List<PedEntry> getAllEntries() throws IOException
	{
		if (filteredEntries != null) return filteredEntries;
		else
		{
			if (entries == null) entries = super.getAllEntries();
			return entries;
		}
	}

	@Override
	public List<PedEntry> getEntries(long from, long to) throws IOException
	{
		return getAllEntries().subList((int) from, (int) to);
	}

	@Override
	public long getNrOfElements() throws IOException
	{
		return getAllEntries().size();
	}

	private List<PedEntry> filter(List<PedEntry> entries, QueryRule filter, List<String> snpNames)
	{
		List<PedEntry> filtered = new ArrayList<PedEntry>();

		for (PedEntry entry : entries)
		{
			// TODO do not hardcode columnnames
			// TODO implement more then just 'equals'
			if (filter.getOperator() == Operator.EQUALS)
			{
				if (filter.getField().equals("IndividualID") && entry.getIndividual().equals(filter.getValue()))
				{
					filtered.add(entry);
				}
				else if (filter.getField().equals("FamilyID") && entry.getFamily().equals(filter.getValue()))
				{
					filtered.add(entry);
				}
				else if (filter.getField().equals("FatherID") && entry.getFather().equals(filter.getValue()))
				{
					filtered.add(entry);
				}
				else if (filter.getField().equals("MotherID") && entry.getMother().equals(filter.getValue()))
				{
					filtered.add(entry);
				}
				else if (filter.getField().equals("Sex")
						&& Byte.valueOf(entry.getSex()).toString().equals(filter.getValue()))
				{
					filtered.add(entry);
				}
				else if (filter.getField().equals("Phenotype")
						&& Double.valueOf(entry.getPhenotype()).toString().equals(filter.getValue()))
				{
					filtered.add(entry);
				}
				else
				{
					boolean found = false;
					int size = snpNames.size();
					for (int i = 0; i < size && !found; i++)
					{
						if (snpNames.get(i).equals(filter.getField())
								&& entry.getBialleles().get(i).toString().equals(filter.getValue()))
						{
							filtered.add(entry);
							found = true;
						}
					}
				}

			}
		}

		return filtered;
	}
}
