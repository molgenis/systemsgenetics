package org.molgenis.vcf.meta;

import java.util.Map;

public class VcfMetaFilter extends VcfMetaEntry
{
	private static final String KEY_ID = "ID";
	private static final String KEY_DESCRIPTION = "Description";

	public VcfMetaFilter(Map<String, String> properties) {
		super(properties);
	}

	@Override
	public String getName()
	{
		return "FILTER";
	}
	
	public String getId()
	{
		return properties.get(KEY_ID);
	}

	public String getDescription()
	{
		return properties.get(KEY_DESCRIPTION);
	}
}
