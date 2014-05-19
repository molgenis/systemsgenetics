package org.molgenis.vcf.meta;

import static org.testng.Assert.assertEquals;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import org.testng.annotations.Test;

public class VcfMetaInfoTest
{

	@Test(expectedExceptions = IllegalArgumentException.class)
	public void VcfMetaInfo()
	{
		new VcfMetaInfo(null);
	}

	@Test
	public void getDescription()
	{
		Map<String, String> properties = new HashMap<String, String>();
		properties.put("Description", "nice description");
		VcfMetaInfo vcfMetaInfo = new VcfMetaInfo(properties);
		assertEquals(vcfMetaInfo.getDescription(), "nice description");
	}

	@Test
	public void getId()
	{
		Map<String, String> properties = new HashMap<String, String>();
		properties.put("ID", "id");
		VcfMetaInfo vcfMetaInfo = new VcfMetaInfo(properties);
		assertEquals(vcfMetaInfo.getId(), "id");
	}

	@Test
	public void getName()
	{
		assertEquals(new VcfMetaInfo(Collections.<String, String> emptyMap()).getName(), "INFO");
	}

	@Test
	public void getNumber()
	{
		Map<String, String> properties = new HashMap<String, String>();
		properties.put("Number", "123");
		VcfMetaInfo vcfMetaInfo = new VcfMetaInfo(properties);
		assertEquals(vcfMetaInfo.getNumber(), "123");
	}

	@Test
	public void getSource()
	{
		Map<String, String> properties = new HashMap<String, String>();
		properties.put("Source", "src");
		VcfMetaInfo vcfMetaInfo = new VcfMetaInfo(properties);
		assertEquals(vcfMetaInfo.getSource(), "src");
	}

	@Test
	public void getType()
	{
		Map<String, String> properties = new HashMap<String, String>();
		properties.put("Type", "Integer");
		VcfMetaInfo vcfMetaInfo = new VcfMetaInfo(properties);
		assertEquals(vcfMetaInfo.getType(), VcfMetaInfo.Type.INTEGER);
	}

	@Test
	public void getVersion()
	{
		Map<String, String> properties = new HashMap<String, String>();
		properties.put("Version", "v1");
		VcfMetaInfo vcfMetaInfo = new VcfMetaInfo(properties);
		assertEquals(vcfMetaInfo.getVersion(), "v1");
	}
}
