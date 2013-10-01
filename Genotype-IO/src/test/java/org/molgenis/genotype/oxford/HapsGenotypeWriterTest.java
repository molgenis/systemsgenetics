package org.molgenis.genotype.oxford;

import org.molgenis.genotype.oxford.HapsGenotypeWriter;
import org.molgenis.genotype.oxford.HapsGenotypeData;
import static org.testng.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;

import org.apache.commons.io.FileUtils;
import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.ResourceTest;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

public class HapsGenotypeWriterTest extends ResourceTest
{
	private GenotypeData genotypeData;
	private HapsGenotypeWriter writer;

	@BeforeClass
	public void setUp() throws IOException, URISyntaxException
	{
		genotypeData = new HapsGenotypeData(getTestImpute2Haps(), getTestImpute2Sample());
		writer = new HapsGenotypeWriter(genotypeData);
	}

	@Test
	public void write() throws IOException, URISyntaxException
	{
		try
		{
			writer.write("write-test");
            
			assertEquals(FileUtils.readFileToString(new File("write-test.haps")),
					FileUtils.readFileToString(getTestImpute2Haps()));

			//TODO fix this test. 
//			assertEquals(FileUtils.readFileToString(new File("write-test.sample")),
//					FileUtils.readFileToString(getTestImpute2Sample()));
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
