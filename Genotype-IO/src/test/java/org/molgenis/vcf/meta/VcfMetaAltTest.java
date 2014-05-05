package org.molgenis.vcf.meta;

import static org.testng.Assert.assertEquals;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import org.testng.annotations.Test;

public class VcfMetaAltTest
{
	@Test(expectedExceptions = IllegalArgumentException.class)
	public void VcfMetaAlt()
	{
		new VcfMetaAlt(null);
	}

	@Test
	public void getDescription()
	{
		Map<String, String> properties = new HashMap<String, String>();
		properties.put("Description", "nice description");
		VcfMetaAlt vcfMetaAlt = new VcfMetaAlt(properties);
		assertEquals(vcfMetaAlt.getDescription(), "nice description");
	}

	@Test
	public void getId()
	{
		Map<String, String> properties = new HashMap<String, String>();
		properties.put("ID", "id");
		VcfMetaAlt vcfMetaAlt = new VcfMetaAlt(properties);
		assertEquals(vcfMetaAlt.getId(), "id");
	}

	@Test
	public void getName()
	{
		assertEquals(new VcfMetaAlt(Collections.<String, String>emptyMap()).getName(), "ALT");
	}
}
