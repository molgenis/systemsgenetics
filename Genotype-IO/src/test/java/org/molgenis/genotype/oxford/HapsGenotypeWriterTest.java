package org.molgenis.genotype.oxford;

import static org.testng.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;

import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.ResourceTest;
import org.molgenis.genotype.util.GenotypeDataCompareTool;
import static org.testng.Assert.assertTrue;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

public class HapsGenotypeWriterTest extends ResourceTest
{
	private GenotypeData genotypeData;
	private HapsGenotypeWriter writer;
	
	private File tmpOutputFolder;
	private String fileSep = System.getProperty("file.separator");

	@BeforeClass
	public void setUp() throws IOException, URISyntaxException
	{
		
		File tmpDir = new File(System.getProperty("java.io.tmpdir"));

		DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
		Date date = new Date();

		tmpOutputFolder = new File(tmpDir, "GenotypeHapsTest_" + dateFormat.format(date));


		Runtime.getRuntime().addShutdownHook(new Thread() {
			@Override
			public void run() {
				System.out.println("Removing tmp dir and files");
				for (File file : tmpOutputFolder.listFiles()) {
					System.out.println(" - Deleting: " + file.getAbsolutePath());
					file.delete();
				}
				System.out.println(" - Deleting: " + tmpOutputFolder.getAbsolutePath());
				tmpOutputFolder.delete();
			}
		});

		tmpOutputFolder.mkdir();

		
		genotypeData = new HapsGenotypeData(getTestImpute2Haps(), getTestImpute2Sample());
		writer = new HapsGenotypeWriter(genotypeData);
	}

	@Test
	public void write() throws IOException, URISyntaxException
	{
	
		writer.write(tmpOutputFolder.getAbsolutePath() + fileSep + "test");

		GenotypeData genotypeDataWritten = RandomAccessGenotypeDataReaderFormats.SHAPEIT2.createGenotypeData(tmpOutputFolder.getAbsolutePath() + fileSep + "test", 0);

		assertTrue(GenotypeDataCompareTool.same(genotypeData, genotypeDataWritten));
		
		assertEquals(genotypeData.getSamples().get(0).getMissingRate(), 0.25, 0.000001);
		assertEquals(genotypeData.getSamples().get(1).getMissingRate(), 0, 0.000001);
		assertEquals(genotypeData.getSamples().get(2).getMissingRate(), 0, 0.000001);
		assertEquals(genotypeData.getSamples().get(3).getMissingRate(), 0, 0.000001);
		
		assertEquals(genotypeData.getSamples().get(0).getAnnotationValues().get("bin1"), Boolean.TRUE);
		assertEquals(genotypeData.getSamples().get(1).getAnnotationValues().get("bin1"), Boolean.FALSE);
		assertEquals(genotypeData.getSamples().get(2).getAnnotationValues().get("bin1"), Boolean.TRUE);
		assertEquals(genotypeData.getSamples().get(3).getAnnotationValues().get("bin1"), Boolean.TRUE);
		
		genotypeData.close();
		
	}
}
