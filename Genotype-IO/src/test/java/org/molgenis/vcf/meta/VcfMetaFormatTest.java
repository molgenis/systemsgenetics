package org.molgenis.vcf.meta;

import static org.testng.Assert.assertEquals;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import org.testng.annotations.Test;

public class VcfMetaFormatTest
{
	@Test(expectedExceptions = IllegalArgumentException.class)
	public void VcfMetaFormat()
	{
		new VcfMetaFormat(null);
	}

	@Test
	public void getDescription()
	{
		Map<String, String> properties = new HashMap<String, String>();
		properties.put("Description", "nice description");
		VcfMetaFormat vcfMetaFormat = new VcfMetaFormat(properties);
		assertEquals(vcfMetaFormat.getDescription(), "nice description");
	}

	@Test
	public void getId()
	{
		Map<String, String> properties = new HashMap<String, String>();
		properties.put("ID", "id");
		VcfMetaFormat vcfMetaFormat = new VcfMetaFormat(properties);
		assertEquals(vcfMetaFormat.getId(), "id");
	}

	@Test
	public void getName()
	{
		assertEquals(new VcfMetaFormat(Collections.<String, String> emptyMap()).getName(), "FORMAT");
	}

	@Test
	public void getNumber()
	{
		Map<String, String> properties = new HashMap<String, String>();
		properties.put("Number", "123");
		VcfMetaFormat vcfMetaFormat = new VcfMetaFormat(properties);
		assertEquals(vcfMetaFormat.getNumber(), "123");
	}

	@Test
	public void getType()
	{
		Map<String, String> properties = new HashMap<String, String>();
		properties.put("Type", "Integer");
		VcfMetaFormat vcfMetaFormat = new VcfMetaFormat(properties);
		assertEquals(vcfMetaFormat.getType(), VcfMetaFormat.Type.INTEGER);
	}
}
