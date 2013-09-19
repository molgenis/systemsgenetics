package org.molgenis.genotype.vcf;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertNotNull;

import org.molgenis.genotype.vcf.VcfHeaderParser;
import org.molgenis.util.tuple.Tuple;
import org.testng.annotations.Test;

public class VcfHeaderParserTest
{
	@Test
	public void testFormatHeader()
	{
		String header = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
		Tuple t = new VcfHeaderParser(header).parse();
		assertNotNull(t);

		assertEquals(t.get("ID"), "GT");
		assertEquals(t.get("Number"), "1");
		assertEquals(t.get("Type"), "String");
		assertEquals(t.get("Description"), "Genotype");
	}

	@Test
	public void testContigHeader()
	{
		String header = "##contig=<ID=20,length=62435964>";
		Tuple t = new VcfHeaderParser(header).parse();

		assertNotNull(t);
		assertEquals(t.getString("ID"), "20");
		assertEquals(t.getLong("length").longValue(), 62435964);
	}
}
