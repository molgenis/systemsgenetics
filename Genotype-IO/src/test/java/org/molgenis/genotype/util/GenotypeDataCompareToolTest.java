package org.molgenis.genotype.util;

import static org.testng.Assert.assertTrue;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URISyntaxException;

import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.ResourceTest;
import org.molgenis.genotype.plink.PedMapGenotypeData;
import org.molgenis.genotype.plink.PedMapGenotypeWriter;
import org.testng.annotations.Test;

public class GenotypeDataCompareToolTest extends ResourceTest
{
	@Test
	public void testPedMap() throws FileNotFoundException, IOException, URISyntaxException
	{
		try
		{
			// Read ped/map
			GenotypeData genotypeData1 = new PedMapGenotypeData(getTestPed(), getTestMap());

			// Write to file
			new PedMapGenotypeWriter(genotypeData1).write("pedmapcomparetest");

			// Read it again
			GenotypeData genotypeData2 = new PedMapGenotypeData("pedmapcomparetest");

			// Should be the same
			assertTrue(GenotypeDataCompareTool.same(genotypeData1, genotypeData2));
		}
		finally
		{
			new File("pedmapcomparetest.ped").delete();
			new File("pedmapcomparetest.map").delete();
		}
	}

}
