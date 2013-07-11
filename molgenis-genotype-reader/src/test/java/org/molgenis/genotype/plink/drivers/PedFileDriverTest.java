package org.molgenis.genotype.plink.drivers;

import static org.testng.AssertJUnit.assertEquals;

import java.io.IOException;

import org.molgenis.genotype.Allele;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

public class PedFileDriverTest extends AbstractResourceTest
{
	private PedFileDriver pedfd;

	@BeforeClass
	public void setup() throws Exception
	{
		System.out.println(getTestResource("/test.ped").getCanonicalPath());
		pedfd = new PedFileDriver(getTestResource("/test.ped"));
	}

	@Test
	public void PED_construct() throws Exception
	{
		assertEquals(9, pedfd.getNrOfElements());
	}

	@Test
	public void PED_getEntries() throws Exception
	{
		// 1 1 0 0 1 1 A A G T
		assertEquals(Allele.C, pedfd.getAllEntries().get(0).getBialleles().get(0).get(0));
		assertEquals(Allele.C, pedfd.getAllEntries().get(0).getBialleles().get(0).get(1));
		assertEquals(Allele.G, pedfd.getAllEntries().get(0).getBialleles().get(1).get(0));
		assertEquals(Allele.G, pedfd.getAllEntries().get(0).getBialleles().get(1).get(1));

		assertEquals(Allele.C, pedfd.getAllEntries().get(1).getBialleles().get(0).get(0));
		assertEquals(Allele.G, pedfd.getAllEntries().get(1).getBialleles().get(0).get(1));
		assertEquals(Allele.A, pedfd.getAllEntries().get(1).getBialleles().get(1).get(0));
		assertEquals(Allele.G, pedfd.getAllEntries().get(1).getBialleles().get(1).get(1));

		assertEquals(Allele.C, pedfd.getAllEntries().get(5).getBialleles().get(0).get(0));
		assertEquals(Allele.C, pedfd.getAllEntries().get(5).getBialleles().get(0).get(1));
		assertEquals(Allele.G, pedfd.getAllEntries().get(5).getBialleles().get(1).get(0));
		assertEquals(Allele.G, pedfd.getAllEntries().get(5).getBialleles().get(1).get(1));

		assertEquals("F1044", pedfd.getAllEntries().get(2).getFamily());
		assertEquals("3", pedfd.getAllEntries().get(2).getFather());
		assertEquals("0", pedfd.getAllEntries().get(2).getMother());
		assertEquals(1.0, pedfd.getAllEntries().get(2).getPhenotype());
		assertEquals("1044", pedfd.getAllEntries().get(2).getIndividual());
		assertEquals(2, pedfd.getAllEntries().get(2).getSex());

		assertEquals("F1044", pedfd.getEntries(2, 3).get(0).getFamily());
		assertEquals(Allele.G, pedfd.getEntries(2, 3).get(0).getBialleles().get(1).get(1));
		assertEquals("F1044", pedfd.getEntries(0, 3).get(2).getFamily());
		assertEquals(Allele.G, pedfd.getEntries(0, 3).get(2).getBialleles().get(1).get(0));
		assertEquals("F1045", pedfd.getEntries(2, 4).get(1).getFamily());
		assertEquals(1.0, pedfd.getEntries(2, 4).get(1).getPhenotype());

		assertEquals(Allele.C, pedfd.getEntries(0, 6).get(4).getBialleles().get(0).get(0));
		assertEquals(Allele.C, pedfd.getEntries(1, 6).get(3).getBialleles().get(0).get(1));
		assertEquals(Allele.G, pedfd.getEntries(2, 6).get(2).getBialleles().get(1).get(0));
		assertEquals(Allele.G, pedfd.getEntries(3, 6).get(1).getBialleles().get(1).get(1));
	}

	@AfterClass
	public void close() throws IOException
	{
		pedfd.close();
	}
}
