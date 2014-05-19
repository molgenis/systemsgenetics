package org.molgenis.vcf.meta;

import static org.testng.Assert.assertEquals;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import org.testng.annotations.Test;

public class VcfMetaSampleTest
{

	@Test(expectedExceptions = IllegalArgumentException.class)
	public void VcfMetaSample()
	{
		new VcfMetaSample(null);
	}

	@Test
	public void getId()
	{
		Map<String, String> properties = new HashMap<String, String>();
		properties.put("ID", "id");
		VcfMetaSample vcfMetaSample = new VcfMetaSample(properties);
		assertEquals(vcfMetaSample.getId(), "id");
	}

	@Test
	public void getName()
	{
		assertEquals(new VcfMetaSample(Collections.<String, String> emptyMap()).getName(), "SAMPLE");
	}
}
