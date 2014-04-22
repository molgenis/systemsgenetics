package org.molgenis.vcf.meta;

import static org.testng.Assert.assertEquals;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import org.testng.annotations.Test;

public class VcfMetaFilterTest {

	@Test(expectedExceptions = IllegalArgumentException.class)
	public void VcfMetaFilter()
	{
		new VcfMetaFilter(null);
	}

	@Test
	public void getDescription()
	{
		Map<String, String> properties = new HashMap<String, String>();
		properties.put("Description", "nice description");
		VcfMetaFilter vcfMetaFilter = new VcfMetaFilter(properties);
		assertEquals(vcfMetaFilter.getDescription(), "nice description");
	}

	@Test
	public void getId()
	{
		Map<String, String> properties = new HashMap<String, String>();
		properties.put("ID", "id");
		VcfMetaFilter vcfMetaFilter = new VcfMetaFilter(properties);
		assertEquals(vcfMetaFilter.getId(), "id");
	}

	@Test
	public void getName()
	{
		assertEquals(new VcfMetaFilter(Collections.<String, String>emptyMap()).getName(), "FILTER");
	}
}
