package org.molgenis.genotype.vcf;

import org.molgenis.util.tuple.Tuple;

public class VcfAlt
{
	private final String id;
	private final String description;

	public VcfAlt(Tuple t)
	{
		id = t.getString("ID");
		description = t.getString("Description");
	}

	public String getId()
	{
		return id;
	}

	public String getDescription()
	{
		return description;
	}

}
