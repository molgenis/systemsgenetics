package org.molgenis.vcf.meta;

import static org.testng.Assert.assertEquals;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import org.testng.annotations.Test;

public class VcfMetaPedigreeTest
{

	@Test(expectedExceptions = IllegalArgumentException.class)
	public void VcfMetaPedigree()
	{
		new VcfMetaPedigree(null);
	}

	@Test
	public void get()
	{
		Map<String, String> properties = new HashMap<String, String>();
		properties.put("name0", "name #0");
		properties.put("name1", "name #1");
		VcfMetaPedigree vcfMetaPedigree = new VcfMetaPedigree(properties);
		assertEquals(vcfMetaPedigree.get("name0"), "name #0");
		assertEquals(vcfMetaPedigree.get("name1"), "name #1");
	}
	
	@Test
	public void getName()
	{
		assertEquals(new VcfMetaPedigree(Collections.<String, String>emptyMap()).getName(), "PEDIGREE");
	}
}
