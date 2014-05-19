package org.molgenis.vcf.meta;

import static org.testng.Assert.assertEquals;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import org.testng.annotations.Test;

public class VcfMetaContigTest
{

	@Test(expectedExceptions = IllegalArgumentException.class)
	public void VcfMetaContig()
	{
		new VcfMetaContig(null);
	}

	@Test
	public void getId()
	{
		Map<String, String> properties = new HashMap<String, String>();
		properties.put("ID", "id");
		VcfMetaContig vcfMetaContig = new VcfMetaContig(properties);
		assertEquals(vcfMetaContig.getId(), "id");
	}

	@Test
	public void getLength()
	{
		Map<String, String> properties = new HashMap<String, String>();
		properties.put("length", "123");
		VcfMetaContig vcfMetaContig = new VcfMetaContig(properties);
		assertEquals(vcfMetaContig.getLength(), Integer.valueOf(123));
	}

	@Test
	public void getName()
	{
		assertEquals(new VcfMetaContig(Collections.<String, String> emptyMap()).getName(), "contig");
	}
}
