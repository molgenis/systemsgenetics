package org.molgenis.genotype.impute2;

import static org.testng.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;

import org.apache.commons.io.FileUtils;
import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.ResourceTest;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

public class Impute2GenotypeWriterTest extends ResourceTest
{
	private GenotypeData genotypeData;
	private Impute2GenotypeWriter writer;

	@BeforeClass
	public void setUp() throws IOException, URISyntaxException
	{
		genotypeData = new Impute2GenotypeData(getTestImpute2Gz(), getTestImpute2GzTbi(), getTestImpute2Sample());
		writer = new Impute2GenotypeWriter(genotypeData);
	}

	@Test
	public void write() throws IOException, URISyntaxException
	{
		try
		{
			writer.write("write-test");

            //TODO fix this test. Did not work for Marc Jan
            
			//assertEquals(FileUtils.readFileToString(new File("write-test.haps")),
			//		FileUtils.readFileToString(getTestImpute2Haps()));

			//assertEquals(FileUtils.readFileToString(new File("write-test.sample")),
			//		FileUtils.readFileToString(getTestImpute2Sample()));
		}
		finally
		{
			File f = new File("write-test.haps");
			if (f.exists())
			{
				f.delete();
			}

			f = new File("write-test.sample");
			if (f.exists())
			{
				f.delete();
			}
		}
	}
}
