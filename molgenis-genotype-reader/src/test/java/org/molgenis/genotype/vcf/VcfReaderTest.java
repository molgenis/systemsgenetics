package org.molgenis.genotype.vcf;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertFalse;
import static org.testng.Assert.assertNotNull;
import static org.testng.Assert.assertNull;
import static org.testng.Assert.assertTrue;

import java.io.FileInputStream;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.Iterator;
import java.util.List;

import org.molgenis.genotype.ResourceTest;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

public class VcfReaderTest extends ResourceTest
{
	private VcfReader reader;

	@BeforeMethod
	public void setUp() throws IOException, URISyntaxException
	{
		reader = new VcfReader(new FileInputStream(getTestVcf1()));
	}

	@AfterMethod
	public void shutDown() throws IOException
	{
		reader.close();
	}

	@Test
	public void testGetSampleNames() throws IOException
	{
		List<String> sampleNames = reader.getSampleNames();

		assertNotNull(sampleNames);
		assertEquals(sampleNames.size(), 1);
		assertTrue(sampleNames.contains("test_S0_L001_R1_001_converted_Unique_Output.pjt"));
	}

	@Test
	public void testVcfIterator()
	{
		int count = 0;
		Iterator<VcfRecord> it = reader.recordIterator();
		while (it.hasNext())
		{
			it.next();
			count++;
		}

		assertEquals(count, 7);
	}

	@Test
	public void testGetInfos()
	{
		List<VcfInfo> infos = reader.getInfos();
		assertNotNull(infos);
		assertEquals(infos.size(), 20);

		// ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
		VcfInfo info = infos.get(0);
		assertEquals(info.getId(), "NS");
		assertEquals(info.getNumber(), "1");
		assertEquals(info.getType(), VcfInfo.Type.INTEGER);
		assertEquals(info.getDescription(), "Number of Samples With Data");
	}

	@Test
	public void testGetFormats()
	{
		List<VcfFormat> formats = reader.getFormats();
		assertNotNull(formats);
		assertEquals(formats.size(), 6);

		// ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
		VcfFormat format = formats.get(0);
		assertEquals(format.getId(), "GT");
		assertEquals(format.getNumber(), "1");
		assertEquals(format.getType(), VcfFormat.Type.STRING);
		assertEquals(format.getDescription(), "Genotype");
	}

	@Test
	public void testGetContigs()
	{
		List<VcfContig> contigs = reader.getContigs();
		assertNotNull(contigs);
		assertEquals(contigs.size(), 1);

		VcfContig contig = contigs.get(0);
		assertEquals(contig.getId(), "1");
		assertEquals(contig.getLength().longValue(), 249240621);
	}

	@Test
	public void testVcfRecord() throws IOException
	{
		Iterator<VcfRecord> it = reader.recordIterator();
		assertNotNull(it);
		assertTrue(it.hasNext());

		VcfRecord record = it.next();
		assertEquals(record.getChrom(), "1");
		assertEquals(record.getPos(), Integer.valueOf(565286));
		assertNotNull(record.getId());
		assertEquals(record.getId().size(), 1);
		assertEquals(record.getId().get(0), "rs1578391");
		assertEquals(record.getRef(), "C");
		assertNotNull(record.getAlt());
		assertEquals(record.getAlt().size(), 1);
		assertEquals(record.getAlt().get(0), "T");
		assertNull(record.getQual());
		assertNotNull(record.getFilter());
		assertEquals(record.getFilter().size(), 1);
		assertEquals(record.getFilter().get(0), "flt");
		assertEquals(record.getSampleValue("test_S0_L001_R1_001_converted_Unique_Output.pjt", "CONFS"),
				"5.300,5.300,1.000,1.000,1.000,1.000,1.000");

		VcfSampleGenotype geno = record.getSampleGenotype("test_S0_L001_R1_001_converted_Unique_Output.pjt");
		assertNotNull(geno);

		assertNotNull(geno.getPhasing());
		assertEquals(geno.getPhasing().size(), 1);
		assertFalse(geno.getPhasing().get(0));

		List<String> sampleVariants = geno.getSamleVariants(record.getAlleles());
		assertNotNull(sampleVariants);
		assertEquals(sampleVariants.size(), 2);
		assertEquals(sampleVariants.get(0), "T");
		assertEquals(sampleVariants.get(1), "T");
	}

	@Test
	public void testGetSamples()
	{
		List<VcfSample> samples = reader.getSamples();
		assertNotNull(samples);
		assertEquals(samples.size(), 1);

		VcfSample sample = samples.get(0);
		assertEquals(sample.getId(), "test_S0_L001_R1_001_converted_Unique_Output.pjt");

		assertNotNull(sample.getDescription());
		assertEquals(sample.getDescription().size(), 1);
		assertEquals(sample.getDescription().get(0), "S1");

		assertNotNull(sample.getGenomes());
		assertEquals(sample.getGenomes().size(), 1);
		assertEquals(sample.getGenomes().get(0), "G1");

		assertNotNull(sample.getMixture());
		assertEquals(sample.getMixture().size(), 1);
		assertEquals(sample.getMixture().get(0), "N1");
	}

	@Test
	public void testGetAlts()
	{
		List<VcfAlt> alts = reader.getAlts();
		assertNotNull(alts);
		assertEquals(alts.size(), 1);
		VcfAlt alt = alts.get(0);
		assertEquals(alt.getId(), "DEL");
		assertEquals(alt.getDescription(), "Deletion");
	}
}
