package org.molgenis.vcf.meta;

import java.util.Map;

public class VcfMetaContig extends VcfMetaEntry
{
	private static final String KEY_ID = "ID";
	private static final String KEY_LENGTH = "length";

	public VcfMetaContig(Map<String, String> properties) {
		super(properties);
	}

	@Override
	public String getName()
	{
		return "contig";
	}
	
	public String getId()
	{
		return properties.get(KEY_ID);
	}

	public Integer getLength()
	{
		return Integer.valueOf(properties.get(KEY_LENGTH));
	}
}
