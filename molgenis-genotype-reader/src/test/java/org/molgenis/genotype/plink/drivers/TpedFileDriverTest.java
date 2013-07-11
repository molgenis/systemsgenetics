package org.molgenis.genotype.plink.drivers;

import static org.testng.AssertJUnit.assertEquals;

import java.io.IOException;
import java.util.List;

import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.plink.datatypes.TpedEntry;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

public class TpedFileDriverTest extends AbstractResourceTest
{
	private TpedFileDriver tpedfd;

	@BeforeClass
	public void setup() throws Exception
	{
		tpedfd = new TpedFileDriver(getTestResource("/test.tped"));
	}

	@Test
	public void TPED_construct() throws Exception
	{
		assertEquals(2, tpedfd.getNrOfElements());
	}

	@Test
	public void TPED_getEntries() throws Exception
	{
		List<TpedEntry> entries = tpedfd.getAllEntries();

		List<Alleles> bialleles = entries.get(0).getBialleles();
		assertEquals(6, bialleles.size());
		assertEquals(Allele.A, bialleles.get(0).get(0));
		assertEquals(Allele.A, bialleles.get(0).get(1));
		assertEquals(Allele.A, bialleles.get(1).get(0));
		assertEquals(Allele.C, bialleles.get(1).get(1));
		assertEquals(Allele.C, bialleles.get(2).get(0));
		assertEquals(Allele.C, bialleles.get(2).get(1));
		assertEquals(Allele.A, bialleles.get(3).get(0));
		assertEquals(Allele.C, bialleles.get(3).get(1));
		assertEquals(Allele.C, bialleles.get(4).get(0));
		assertEquals(Allele.C, bialleles.get(4).get(1));
		assertEquals(Allele.C, bialleles.get(5).get(0));
		assertEquals(Allele.C, bialleles.get(5).get(1));
	}

	@AfterClass
	public void close() throws IOException
	{
		tpedfd.close();
	}
}
