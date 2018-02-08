/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.bgen;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import org.molgenis.genotype.ResourceTest;
import static org.testng.Assert.*;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

/**
 *
 * @author Patrick Deelen
 */
public class BgenGenotypeDataTest extends ResourceTest {

	private BgenGenotypeData genotypeData;

//	@BeforeClass
//	public void beforeClass() throws IOException, URISyntaxException {
//		genotypeData = new BgenGenotypeData(getTest3Bgen(), getTest3Sample());
//	}
	@Test
	public void testSomeMethod() {
		// TODO review the generated test code and remove the default call to fail.
	}

	@Test
	public void testReader_1_1() throws URISyntaxException, IOException {
		genotypeData = new BgenGenotypeData(getTest3Bgen(), getTest3Sample());
	}

	@Test
	public void testReader_1_2() throws URISyntaxException, IOException {
		//File bgenixFile = getTestResourceFile("/bgenExamples/complex.bgen");
		//genotypeData = new BgenGenotypeData(bgenixFile, null);
		File bgenixFile = getTestResourceFile("/bgenExamples/example.16bits.bgen");
		genotypeData = new BgenGenotypeData(bgenixFile, null);
	}

	@Test
	public void testReader_1_3() throws URISyntaxException, IOException {
		File bgenixFile = getTestResourceFile("/bgenExamples/example.16bits.zstd.bgen");
		genotypeData = new BgenGenotypeData(bgenixFile, null);
	}
}
