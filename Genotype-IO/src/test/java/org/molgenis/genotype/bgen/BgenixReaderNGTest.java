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
public class BgenixReaderNGTest extends ResourceTest{
	
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

		ResultSet result = instance.getVariantsChromosome(chr);
		
		System.out.println("Test");
		
		while(result.next()){
			System.out.print(result.getString("chromosome"));
			System.out.print("\t");
			System.out.print(result.getString("position"));
			System.out.print("\t");
			System.out.print(result.getString("rsid"));
			System.out.print("\t");
			System.out.print(result.getString("allele1"));
			System.out.print("\t");
			System.out.print(result.getString("allele2"));
			System.out.println();
		}
		
		//assertEquals(result, expResult);
		// TODO review the generated test code and remove the default call to fail.
		//fail("The test case is a prototype.");
	}
	
}
