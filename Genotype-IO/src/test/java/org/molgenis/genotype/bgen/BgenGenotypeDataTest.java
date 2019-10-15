/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.bgen;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.junit.Rule;
import org.junit.rules.TemporaryFolder;
import org.molgenis.genotype.ResourceTest;
import org.testng.annotations.Test;

/**
 *
 * @author Patrick Deelen
 */
public class BgenGenotypeDataTest extends ResourceTest {

	private BgenGenotypeData genotypeData;

	@Rule
	public TemporaryFolder folder = new TemporaryFolder();

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
		//genotypeData = new BgenGenotypeData(getTest3Bgen(), getTest3Sample());
	}

	@Test
	public void testReader_1_2() throws URISyntaxException, IOException {
		//File bgenFile = getTestResourceFile("/bgenExamples/complex.bgen");
		//genotypeData = new BgenGenotypeData(bgenFile, null);
		folder.create();
		File bgenFile = getTestResourceFile("/bgenExamples/example.25bits.bgen");
		Path target = Paths.get(folder.getRoot().toString(), bgenFile.getName());
		Files.copy(bgenFile.toPath(), target);
		genotypeData = new BgenGenotypeData(target.toFile(), null);
	}

	@Test
	public void testReader_1_3() throws URISyntaxException, IOException {
		folder.create();
		File bgenFile = getTestResourceFile("/bgenExamples/example.16bits.zstd.bgen");
		Path target = Paths.get(folder.getRoot().toString(), bgenFile.getName());
		Files.copy(bgenFile.toPath(), target);
		genotypeData = new BgenGenotypeData(target.toFile(), null);
	}
}
