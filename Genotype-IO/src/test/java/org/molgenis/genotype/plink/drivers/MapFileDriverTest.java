package org.molgenis.genotype.plink.drivers;

import static org.testng.AssertJUnit.assertEquals;

import java.io.IOException;

import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

public class MapFileDriverTest extends AbstractResourceTest
{
	private MapFileDriver mapfd;

	@BeforeClass
	public void setup() throws Exception
	{
		mapfd = new MapFileDriver(getTestResource("/test.map"));
	}

	@Test
	public void MAP_construct() throws Exception
	{
		assertEquals(10, mapfd.getNrOfElements());
	}

	@Test
	public void MAP_getEntries() throws Exception
	{
		assertEquals(1, mapfd.getEntries(0, 1).size());
		assertEquals(1, mapfd.getEntries(1, 2).size());
		assertEquals(2, mapfd.getEntries(0, 2).size());
		assertEquals(10, mapfd.getAllEntries().size());
		assertEquals("rs11089130", mapfd.getEntries(0, 1).get(0).getSNP());
		assertEquals(0.0, mapfd.getEntries(1, 2).get(0).getcM());
		assertEquals("22", mapfd.getEntries(1, 2).get(0).getChromosome());
		assertEquals(14431347, mapfd.getAllEntries().get(0).getBpPos());
		assertEquals(14432618, mapfd.getAllEntries().get(1).getBpPos());
		assertEquals("rs738829", mapfd.getAllEntries().get(1).getSNP());
	}

	@AfterClass
	public void close() throws IOException
	{
		mapfd.close();
	}
}
