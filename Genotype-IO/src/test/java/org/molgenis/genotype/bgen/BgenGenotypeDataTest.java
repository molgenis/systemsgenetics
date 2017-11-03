/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.bgen;

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

	@BeforeClass
	public void beforeClass() throws IOException, URISyntaxException {
		genotypeData = new BgenGenotypeData(getTest3Bgen(), getTest3Sample());
	}
	
	@Test
	public void testSomeMethod() {
		// TODO review the generated test code and remove the default call to fail.
	
	}
}