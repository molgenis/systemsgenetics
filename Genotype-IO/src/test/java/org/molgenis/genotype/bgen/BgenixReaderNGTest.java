/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.bgen;

import java.io.File;
import java.net.URISyntaxException;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Arrays;
import static org.testng.Assert.*;
import org.testng.annotations.AfterClass;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;
import org.molgenis.genotype.ResourceTest;

/**
 *
 * @author patri
 */
public class BgenixReaderNGTest extends ResourceTest {

	public BgenixReaderNGTest() {
	}

	@BeforeClass
	public static void setUpClass() throws Exception {
	}

	@AfterClass
	public static void tearDownClass() throws Exception {
	}

	@BeforeMethod
	public void setUpMethod() throws Exception {
	}

	@AfterMethod
	public void tearDownMethod() throws Exception {
	}

	/**
	 * Test of getVariantsChromosome method, of class BgenixReader.
	 */
	@Test
	public void testGetVariantsChromosome() throws URISyntaxException, SQLException {
		File bgenixFile = getTestResourceFile("/bgenExamples/complex.bgen.bgi");

		String chr = "01";
		BgenixReader instance = new BgenixReader(bgenixFile);

//		ResultSet result = instance.getVariantsChromosome(chr);
//
//
//		while (result.next()) {
//			System.out.print(result.getString("chromosome"));
//			System.out.print("\t");
//			System.out.print(result.getString("position"));
//			System.out.print("\t");
//			System.out.print(result.getString("rsid"));
//			System.out.print("\t");
//			System.out.print(result.getString("number_of_alleles"));
//			System.out.print("\t");
//			System.out.print(result.getString("allele1"));
//			System.out.print("\t");
//			System.out.print(result.getString("allele2"));
//			System.out.print("\t");
//			System.out.print(result.getString("file_start_position"));
//			System.out.print("\t");
//			System.out.print(result.getString("size_in_bytes"));
//			System.out.println();
//		}

		//assertEquals(result, expResult);
		// TODO review the generated test code and remove the default call to fail.
		//fail("The test case is a prototype.");
	}

	/**
	 * Test of createNewConnection method, of class BgenixReader.
	 */
	@Test
	public void testCreateNewConnection() {

	}

	/**
	 * Test of getVariantsPostion method, of class BgenixReader.
	 */
	@Test
	public void testGetVariantsPostion() throws URISyntaxException, SQLException {
		File bgenixFile = getTestResourceFile("/bgenExamples/complex.bgen.bgi");

		String chr = "01";
		int pos = 2;
		BgenixReader instance = new BgenixReader(bgenixFile);

//		ResultSet result = instance.getVariantsPostion(chr, pos);
//
//		while (result.next()) {
//			System.out.print(result.getString("chromosome"));
//			System.out.print("\t");
//			System.out.print(result.getString("position"));
//			System.out.print("\t");
//			System.out.print(result.getString("rsid"));
//			System.out.print("\t");
//			System.out.print(result.getString("number_of_alleles"));
//			System.out.print("\t");
//			System.out.print(result.getString("allele1"));
//			System.out.print("\t");
//			System.out.print(result.getString("allele2"));
//			System.out.print("\t");
//			System.out.print(result.getString("file_start_position"));
//			System.out.print("\t");
//			System.out.print(result.getString("size_in_bytes"));
//			System.out.println();
//		}
	}

	/**
	 * Test of getVariantsRange method, of class BgenixReader.
	 */
	@Test
	public void testGetVariantsRange() throws URISyntaxException, SQLException {
		File bgenixFile = getTestResourceFile("/bgenExamples/complex.bgen.bgi");

		String chr = "01";
		int from = 2;
		int to = 4;
		BgenixReader instance = new BgenixReader(bgenixFile);

//		ResultSet result = instance.getVariantsRange(chr, from, to);
//
//		while (result.next()) {
//			System.out.print(result.getString("chromosome"));
//			System.out.print("\t");
//			System.out.print(result.getString("position"));
//			System.out.print("\t");
//			System.out.print(result.getString("rsid"));
//			System.out.print("\t");
//			System.out.print(result.getString("number_of_alleles"));
//			System.out.print("\t");
//			System.out.print(result.getString("allele1"));
//			System.out.print("\t");
//			System.out.print(result.getString("allele2"));
//			System.out.print("\t");
//			System.out.print(result.getString("file_start_position"));
//			System.out.print("\t");
//			System.out.print(result.getString("size_in_bytes"));
//			System.out.println();
//		}
	}

	/**
	 * Test of getMetadata method, of class BgenixReader.
	 */
	@Test
	public void testGetMetadata() throws Exception {
		File bgenixFile = getTestResourceFile("/bgenExamples/complex.bgen.bgi");
		
		BgenixReader bgenixReader = new BgenixReader(bgenixFile);
		
		BgenixMetadata metaData = bgenixReader.getMetadata();
		
		assertEquals(metaData.getFileName(), "complex.bgen");
		assertEquals(metaData.getFileSize(), 835);
		assertEquals(metaData.getIndexCreationTime(), 2017);
		assertEquals(metaData.getLastWriteTime(), 1464771094);
				
	}

}
