package org.molgenis.genotype.plink.drivers;

import static org.testng.AssertJUnit.assertEquals;

import java.io.IOException;

import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

/**
 * This test looks like the PED driver test except there is no allele information.
 * @author jvelde
 *
 */
public class FamFileDriverTest extends AbstractResourceTest
{
	private FamFileDriver famfd;

	@BeforeClass
	public void setup() throws Exception
	{
		famfd = new FamFileDriver(getTestResource("/test.fam"));
	}

	@Test
	public void FAM_construct() throws Exception
	{
		assertEquals(9, famfd.getNrOfElements());
	}

	@Test
	public void FAM_getEntries() throws Exception
	{		
		assertEquals("F1044", famfd.getAllEntries().get(2).getFamily());
		assertEquals("3", famfd.getAllEntries().get(2).getFather());
		assertEquals("0", famfd.getAllEntries().get(2).getMother());
		assertEquals(1.0, famfd.getAllEntries().get(2).getPhenotype());
		assertEquals("1044", famfd.getAllEntries().get(2).getIndividual());
		assertEquals(2, famfd.getAllEntries().get(2).getSex());
		assertEquals("F1044", famfd.getEntries(2, 3).get(0).getFamily());
		assertEquals("F1044", famfd.getEntries(0, 3).get(2).getFamily());
		assertEquals("F1045", famfd.getEntries(2, 4).get(1).getFamily());
		assertEquals(1.0, famfd.getEntries(2, 4).get(1).getPhenotype());
	}

	@AfterClass
	public void close() throws IOException
	{
		famfd.close();
	}
}
