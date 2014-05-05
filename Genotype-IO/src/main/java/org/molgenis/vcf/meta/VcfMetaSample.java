package org.molgenis.vcf.meta;

import java.util.Map;

public class VcfMetaSample extends VcfMetaEntry
{
	private static final String KEY_ID = "ID";
	
	public VcfMetaSample(Map<String, String> properties) {
		super(properties);
	}

	@Override
	public String getName()
	{
		return "SAMPLE";
	}
	
	public String getId()
	{
		return properties.get(KEY_ID);
	}
}
