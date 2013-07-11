package org.molgenis.genotype.impute2;

import static org.testng.Assert.assertEquals;

import java.util.Arrays;

import org.testng.annotations.Test;

public class HapsLineParserTest
{

	@Test
	public void parse()
	{
		String line = "7	SNP1	123	A	G	0	0	1	0	0	0	1	1";
		HapsEntry entry = HapsLineParser.parse(line);
		assertEquals(entry.getChromosomeNumber(), 7);
		assertEquals(entry.getFirstAllele(), "A");
		assertEquals(entry.getPosition(), 123);
		assertEquals(entry.getSecondAllele(), "G");
		assertEquals(entry.getSnpId(), "SNP1");
		assertEquals(entry.getSampleAlleles(), Arrays.asList(new String[]
		{ "0", "0" }, new String[]
		{ "1", "0" }, new String[]
		{ "0", "0" }, new String[]
		{ "1", "1" }));
	}
}
