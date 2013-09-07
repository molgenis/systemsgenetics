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
import static org.testng.Assert.*;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

/**
 *
 * @author Patrick Deelen
 */
public class PedMapGenotypeWriterNGTest extends ResourceTest {

	private File tmpOutputFolder;
	private String fileSep = System.getProperty("file.separator");

	public PedMapGenotypeWriterNGTest() {

		File tmpDir = new File(System.getProperty("java.io.tmpdir"));

		DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
		Date date = new Date();

		tmpOutputFolder = new File(tmpDir, "GenotypePedMapWriterTest_" + dateFormat.format(date));


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


		System.out.println("Temp folder with output of this test: " + tmpOutputFolder.getAbsolutePath());

	}

	@BeforeMethod
	public void setUpMethod() throws Exception {
	}

	/**
	 * Test of write method, of class PedMapGenotypeWriter.
	 */
	@Test
	public void testWrite() throws Exception {
	
		GenotypeData genotypeData = new PedMapGenotypeData(getTestPed(), getTestMap());

		PedMapGenotypeWriter writer = new PedMapGenotypeWriter(genotypeData);

		writer.write(tmpOutputFolder.getAbsolutePath() + fileSep + "test");

		GenotypeData genotypeDataWritten = RandomAccessGenotypeDataReaderFormats.PED_MAP.createGenotypeData(tmpOutputFolder.getAbsolutePath() + fileSep + "test", 0);

		assertTrue(GenotypeDataCompareTool.same(genotypeData, genotypeDataWritten));
		
	}

}
