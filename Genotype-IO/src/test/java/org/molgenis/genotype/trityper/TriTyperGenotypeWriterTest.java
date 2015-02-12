/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.trityper;

import java.io.File;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.GenotypeWriter;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.ResourceTest;
import org.molgenis.genotype.oxford.GenGenotypeData;
import org.molgenis.genotype.sampleFilter.SampleIdIncludeFilter;
import org.molgenis.genotype.util.GenotypeDataCompareTool;
import static org.testng.Assert.*;

import org.molgenis.genotype.variant.GeneticVariant;
import org.testng.annotations.Test;

/**
 *
 * @author Patrick Deelen
 */
public class TriTyperGenotypeWriterTest extends ResourceTest{
	
	private File tmpOutputFolder;
	private String fileSep = System.getProperty("file.separator");
	
	public TriTyperGenotypeWriterTest() {
		
		File tmpDir = new File(System.getProperty("java.io.tmpdir"));

		DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
		Date date = new Date();

		tmpOutputFolder = new File(tmpDir, "GenotypeTriTyperWriterTest_" + dateFormat.format(date));


		Runtime.getRuntime().addShutdownHook(new Thread() {
			@Override
			public void run() {
				System.out.println("Removing tmp dir and files");
				for (File file : tmpOutputFolder.listFiles()) {
					System.out.println(" - Deleting: " + file.getAbsolutePath());
                    //Why?
				//	file.delete();
				}
				System.out.println(" - Deleting: " + tmpOutputFolder.getAbsolutePath());
//				tmpOutputFolder.delete();
			}
		});

		tmpOutputFolder.mkdir();

		System.out.println("Temp folder with output of this test: " + tmpOutputFolder.getAbsolutePath());
		
	}

	@Test
	public void writeTriTyper() throws Exception{
		
		GenotypeData original = new TriTyperGenotypeData(getTriTyperFolder().getAbsolutePath());
		
		GenotypeWriter copyWriter = new TriTyperGenotypeWriter(original);
		
		copyWriter.write(tmpOutputFolder.getAbsolutePath());
		
		GenotypeData copy = new TriTyperGenotypeData(tmpOutputFolder);
		
		assertTrue(GenotypeDataCompareTool.same(original, copy));
		
	}

	@Test
	public void writeTriTyperDosage() throws Exception{

		RandomAccessGenotypeData original = new GenGenotypeData(getTest2Gen(), getTest2Sample());

		GenotypeWriter writer = new TriTyperGenotypeWriter(original);

		writer.write(tmpOutputFolder.getAbsolutePath());

		RandomAccessGenotypeData trityper = new TriTyperGenotypeData(tmpOutputFolder);

		int numSamples = trityper.getSamples().size();

		GeneticVariant variant = original.getSnpVariantByPos("1", 1);
		GeneticVariant originalVariant = original.getSnpVariantByPos("1", 1);

		float[] dosageValues = variant.getSampleDosages();
		float[] originalDosageValues = originalVariant.getSampleDosages();

		System.out.print("\n" + variant.getSequenceName() + ":" + variant.getStartPos());
		for(int i = 0; i < numSamples; i++){
			assertEquals(dosageValues[i], originalDosageValues[i], 0.01);
		}


		variant = original.getSnpVariantByPos("22", 14432918);
		originalVariant = original.getSnpVariantByPos("22", 14432918);

		dosageValues = variant.getSampleDosages();
		originalDosageValues = originalVariant.getSampleDosages();

		System.out.print("\n" + variant.getSequenceName() + ":" + variant.getStartPos());
		for(int i = 0; i < numSamples; i++){
			assertEquals(dosageValues[i], originalDosageValues[i], 0.01);
		}
		
		RandomAccessGenotypeData trityper2 = new TriTyperGenotypeData(tmpOutputFolder, 0, null, new SampleIdIncludeFilter("1042", "1047"));
		
		assertEquals(trityper2.getSamples().size(),2);
		
		variant = trityper2.getSnpVariantByPos("1", 1);
		originalVariant = original.getSnpVariantByPos("1", 1);
		
		float[] dosages = variant.getSampleDosages();
		originalDosageValues = originalVariant.getSampleDosages();
		
		assertEquals(dosages[0], originalDosageValues[0], 0.00001);
		assertEquals(dosages[1], originalDosageValues[5], 0.00001);
        
	}
    
    @Test
    public void writeTriTyperDosage2() throws Exception{

		RandomAccessGenotypeData original = new GenGenotypeData(getTest10Gen(), getTest2Sample());

		GenotypeWriter writer = new TriTyperGenotypeWriter(original);

		writer.write(tmpOutputFolder.getAbsolutePath());

		RandomAccessGenotypeData trityper = new TriTyperGenotypeData(tmpOutputFolder);

		int numSamples = trityper.getSamples().size();

		GeneticVariant variant = original.getSnpVariantByPos("1", 1);
		GeneticVariant originalVariant = original.getSnpVariantByPos("1", 1);

		float[] dosageValues = variant.getSampleDosages();
		float[] originalDosageValues = originalVariant.getSampleDosages();

		System.out.print("\n" + variant.getSequenceName() + ":" + variant.getStartPos());
		for(int i = 0; i < numSamples; i++){
			assertEquals(dosageValues[i], originalDosageValues[i], 0.01);
		}

        //Check allels
//        System.out.print("\n" + variant.getSequenceName() + ":" + variant.getStartPos());
//		for(int i = 0; i < numSamples; i++){
//			assertEquals(dosageValues[i], originalDosageValues[i], 0.01);
//		}
        

		variant = original.getSnpVariantByPos("22", 14432918);
		originalVariant = original.getSnpVariantByPos("22", 14432918);

		dosageValues = variant.getSampleDosages();
		originalDosageValues = originalVariant.getSampleDosages();

		System.out.print("\n" + variant.getSequenceName() + ":" + variant.getStartPos());
		for(int i = 0; i < numSamples; i++){
			assertEquals(dosageValues[i], originalDosageValues[i], 0.01);
		}
		
        //Check allels
//        System.out.print("\n" + variant.getSequenceName() + ":" + variant.getStartPos());
//		for(int i = 0; i < numSamples; i++){
//			assertEquals(dosageValues[i], originalDosageValues[i], 0.01);
//		}
	}
}