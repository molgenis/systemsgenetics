package org.molgenis.vcf.meta;

import java.util.Map;

public class VcfMetaPedigree extends VcfMetaEntry
{
	public VcfMetaPedigree(Map<String, String> properties) {
		super(properties);
	}

	@Override
	public String getName()
	{
		return "PEDIGREE";
	}
}
