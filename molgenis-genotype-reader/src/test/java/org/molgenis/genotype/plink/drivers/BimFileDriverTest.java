package org.molgenis.genotype.plink.drivers;

import static org.testng.AssertJUnit.assertEquals;

import java.io.IOException;

import org.molgenis.genotype.Allele;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

/**
 * Looks like the test for MAP file, but with 2 additional columns for the
 * alleles (encoded in the BED file)
 * 
 * @author jvelde
 * 
 */
public class BimFileDriverTest extends AbstractResourceTest
{
	private BimFileDriver bimfd;

	@BeforeClass
	public void setup() throws Exception
	{
		bimfd = new BimFileDriver(getTestResource("/test.bim"), '\t');
	}

	@Test
	public void BIM_construct() throws Exception
	{
		assertEquals(10, bimfd.getNrOfElements());
	}

	@Test
	public void BIM_getEntries() throws Exception
	{
		assertEquals(1, bimfd.getEntries(0, 1).size());
		assertEquals(1, bimfd.getEntries(1, 2).size());
		assertEquals(2, bimfd.getEntries(0, 2).size());
		assertEquals(10, bimfd.getAllEntries().size());
		assertEquals("rs11089130", bimfd.getEntries(0, 1).get(0).getSNP());
		assertEquals(0.0, bimfd.getEntries(1, 2).get(0).getcM());
		assertEquals("22", bimfd.getEntries(1, 2).get(0).getChromosome());
		assertEquals(14431347, bimfd.getAllEntries().get(0).getBpPos());
		assertEquals(14432618, bimfd.getAllEntries().get(1).getBpPos());
		assertEquals("rs738829", bimfd.getAllEntries().get(1).getSNP());

		// alleles
		assertEquals(Allele.G, bimfd.getAllEntries().get(0).getBiallele().get(0));
		assertEquals(Allele.C, bimfd.getAllEntries().get(0).getBiallele().get(1));
		assertEquals(Allele.A, bimfd.getAllEntries().get(1).getBiallele().get(0));
		assertEquals(Allele.G, bimfd.getAllEntries().get(1).getBiallele().get(1));
		assertEquals(Allele.A, bimfd.getAllEntries().get(5).getBiallele().get(0));
		assertEquals(Allele.C, bimfd.getAllEntries().get(5).getBiallele().get(1));
	}

	@AfterClass
	public void close() throws IOException
	{
		bimfd.close();
	}
}
