/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.bgen;

import java.io.File;
import java.io.RandomAccessFile;
import org.molgenis.genotype.ResourceTest;
import static org.testng.Assert.*;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

/**
 *
 * @author patri
 */
public class BgenixReaderTest extends ResourceTest {

	public BgenixReaderTest() {
	}

	@BeforeClass
	public static void setUpClass() throws Exception {
	}

	@BeforeMethod
	public void setUpMethod() throws Exception {
	}

	@Test
	public void bgenixReaderTest() throws Exception {

		File bgenFile = getTestResourceFile("/bgenExamples/example.16bits.bgen");
		File bgenixFile = getTestResourceFile("/bgenExamples/example.16bits.bgen.bgi");

		BgenixReader bgenix = new BgenixReader(bgenixFile);

		BgenixMetadata bgenixMeta = bgenix.getMetadata();

		System.out.println(bgenixMeta.getFileName());
		System.out.println(bgenixMeta.getFileSize());
		System.out.println(bgenixMeta.getIndexCreationTime());
		System.out.println(bgenixMeta.getLastWriteTime());

		assertEquals(bgenixMeta.getFileName(), "example.16bits.bgen");
		assertEquals(bgenixMeta.getFileSize(), bgenFile.length());

		byte[] snpInfoBuffer = new byte[1000];
		RandomAccessFile bgenReader = new RandomAccessFile(bgenFile, "r");
		bgenReader.read(snpInfoBuffer, 0, snpInfoBuffer.length);
		assertEquals(bgenixMeta.getFirst1000bytes(), snpInfoBuffer);

		//We can't test this since this is flexible.
//		assertEquals(bgenixMeta.getIndexCreationTime(), 2017);
		assertEquals(bgenixMeta.getLastWriteTime(), 1568284398);

		System.out.println(bgenix.getVariantCount());

		BgenixVariantQueryResult variantsChromosome = bgenix.getVariantsChromosome("02");
		assertEquals(variantsChromosome.hasNext(), false);

		variantsChromosome = bgenix.getVariantsChromosome("01");
		assertEquals(variantsChromosome.hasNext(), true);

		variantsChromosome.next();
		variantsChromosome.next();
		BgenixVariantData variant = variantsChromosome.next();

		assertEquals(variant.getChromosome(), "01");
		assertEquals(variant.getRsid(), "RSID_102");
		assertEquals(variant.getPosition(), 2001);
		assertEquals(variant.getNumber_of_alleles(), 2);
		assertEquals(variant.getAllele1(), "A");
		assertEquals(variant.getAllele2(), "G");

		BgenixVariantQueryResult variantsRange = bgenix.getVariantsRange("01", 31000, 32000);
		variant = variantsRange.next();
		assertEquals(variant.getChromosome(), "01");
		assertEquals(variant.getRsid(), "RSID_31");
		assertEquals(variant.getPosition(), 31000);
		assertEquals(variant.getNumber_of_alleles(), 2);
		assertEquals(variant.getAllele1(), "A");
		assertEquals(variant.getAllele2(), "G");

		variant = variantsRange.next();
		assertEquals(variant.getChromosome(), "01");
		assertEquals(variant.getRsid(), "RSID_131");
		assertEquals(variant.getPosition(), 31001);
		assertEquals(variant.getNumber_of_alleles(), 2);
		assertEquals(variant.getAllele1(), "A");
		assertEquals(variant.getAllele2(), "G");
		
		variant = variantsRange.next();
		assertEquals(variant.getChromosome(), "01");
		assertEquals(variant.getRsid(), "RSID_32");
		assertEquals(variant.getPosition(), 32000);
		assertEquals(variant.getNumber_of_alleles(), 2);
		assertEquals(variant.getAllele1(), "A");
		assertEquals(variant.getAllele2(), "G");
		
		assertEquals(variantsRange.hasNext(), false);

		BgenixVariantQueryResult variantsPostion = bgenix.getVariantsPostion("01", 92001);
		variant = variantsPostion.next();
		assertEquals(variant.getChromosome(), "01");
		assertEquals(variant.getRsid(), "RSID_192");
		assertEquals(variant.getPosition(), 92001);
		assertEquals(variant.getNumber_of_alleles(), 2);
		assertEquals(variant.getAllele1(), "A");
		assertEquals(variant.getAllele2(), "G");
		
		assertEquals(variantsPostion.hasNext(), false);
		
	}

}
