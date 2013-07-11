package org.molgenis.genotype.util;

import static org.testng.Assert.assertEquals;

import java.util.Arrays;
import java.util.List;

import org.testng.annotations.Test;

public class UtilsTest
{

	@Test
	public void iteratorToList()
	{
		List<String> list = Arrays.asList("1", "2", null);
		assertEquals(Utils.iteratorToList(list.iterator()), list);
	}

	@Test
	public void swapSnpStrand()
	{
		char[] snps = new char[]
		{ 'A', 'G', 'C', 'T', '0' };
		assertEquals(Utils.swapSnpStrand(snps), new char[]
		{ 'T', 'C', 'G', 'A', '0' });
	}
}
