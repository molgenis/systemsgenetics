/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.plink;

import java.io.File;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.ResourceTest;
import org.molgenis.genotype.util.GenotypeDataCompareTool;
import org.testng.annotations.AfterTest;
import org.testng.annotations.Test;
import static org.testng.Assert.assertTrue;

/**
 *
 * @author Patrick Deelen
 */
public class BedBimFamGenotypeWriterTest extends ResourceTest {

	private File tmpOutputFolder;
	private String fileSep = System.getProperty("file.separator");

	public BedBimFamGenotypeWriterTest() {

		File tmpDir = new File(System.getProperty("java.io.tmpdir"));

		DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
		Date date = new Date();

		tmpOutputFolder = new File(tmpDir, "GenotypeBedBimFamWriterTest_" + dateFormat.format(date));


		Runtime.getRuntime().addShutdownHook(new Thread() {
			@Override
			public void run() {
				System.out.println("Removing tmp dir and files");
				for (File file : tmpOutputFolder.listFiles()) {
					System.out.println(" - Deleting: " + file.getAbsolutePath());
					//file.delete();
				}
				System.out.println(" - Deleting: " + tmpOutputFolder.getAbsolutePath());
				tmpOutputFolder.delete();
			}
		});

		tmpOutputFolder.mkdir();


		System.out.println("Temp folder with output of this test: " + tmpOutputFolder.getAbsolutePath());

	}

	/**
	 * Test of write method, of class BedBimFamGenotypeWriter.
	 */
	@Test
	public void testWrite_String() throws Exception {
		GenotypeData genotypeData = new BedBimFamGenotypeData(getTestBed6(), getTestBim6(), getTestFam6(), 0);

		BedBimFamGenotypeWriter writer = new BedBimFamGenotypeWriter(genotypeData);

		writer.write(tmpOutputFolder.getAbsolutePath() + fileSep + "test6samples");

		GenotypeData genotypeDataWritten = RandomAccessGenotypeDataReaderFormats.PLINK_BED.createGenotypeData(tmpOutputFolder.getAbsolutePath() + fileSep + "test6samples", 0);

		assertTrue(GenotypeDataCompareTool.same(genotypeData, genotypeDataWritten));

	}
	
		/**
	 * Test of write method, of class BedBimFamGenotypeWriter.
	 */
	@Test
	public void testWrite_String2() throws Exception {
		GenotypeData genotypeData = new BedBimFamGenotypeData(getTestBed7(), getTestBim7(), getTestFam7(), 0);

		BedBimFamGenotypeWriter writer = new BedBimFamGenotypeWriter(genotypeData);

		writer.write(tmpOutputFolder.getAbsolutePath() + fileSep + "test7samples");

		GenotypeData genotypeDataWritten = RandomAccessGenotypeDataReaderFormats.PLINK_BED.createGenotypeData(tmpOutputFolder.getAbsolutePath() + fileSep + "test7samples", 0);
		
		assertTrue(GenotypeDataCompareTool.same(genotypeData, genotypeDataWritten));

	}
	
		/**
	 * Test of write method, of class BedBimFamGenotypeWriter.
	 */
	@Test
	public void testWrite_String3() throws Exception {
		GenotypeData genotypeData = new BedBimFamGenotypeData(getTestBed8(), getTestBim8(), getTestFam8(), 0);

		BedBimFamGenotypeWriter writer = new BedBimFamGenotypeWriter(genotypeData);

		writer.write(tmpOutputFolder.getAbsolutePath() + fileSep + "test8samples");

		GenotypeData genotypeDataWritten = RandomAccessGenotypeDataReaderFormats.PLINK_BED.createGenotypeData(tmpOutputFolder.getAbsolutePath() + fileSep + "test8samples", 0);

		assertTrue(GenotypeDataCompareTool.same(genotypeData, genotypeDataWritten));

	}
	
		/**
	 * Test of write method, of class BedBimFamGenotypeWriter.
	 */
	@Test
	public void testWrite_String4() throws Exception {
		GenotypeData genotypeData = new BedBimFamGenotypeData(getTestBed9(), getTestBim9(), getTestFam9(), 0);

		BedBimFamGenotypeWriter writer = new BedBimFamGenotypeWriter(genotypeData);

		writer.write(tmpOutputFolder.getAbsolutePath() + fileSep + "test9samples");

		GenotypeData genotypeDataWritten = RandomAccessGenotypeDataReaderFormats.PLINK_BED.createGenotypeData(tmpOutputFolder.getAbsolutePath() + fileSep + "test9samples", 0);

		assertTrue(GenotypeDataCompareTool.same(genotypeData, genotypeDataWritten));
		
		genotypeData.close();
		genotypeDataWritten.close();

	}

	@AfterTest
	public void removeTempFiles() {
		
		System.out.println("Removing tmp dir and files");
		System.out.println(" - Deleting: " + tmpOutputFolder.getAbsolutePath());
		tmpOutputFolder.deleteOnExit();
		for (File file : tmpOutputFolder.listFiles()) {
			System.out.println(" - Deleting: " + file.getAbsolutePath());
			//file.deleteOnExit();
		}
	}
}
