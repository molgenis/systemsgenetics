/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype;

import java.net.URISyntaxException;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;
import static org.testng.Assert.assertEquals;

/**
 *
 * @author Patrick Deelen
 */
public class RandomAccessGenotypeDataReaderFormatsTest extends ResourceTest{
	
	public RandomAccessGenotypeDataReaderFormatsTest() {
	}

	@BeforeMethod
	public void setUpMethod() throws Exception {
	}


	/**
	 * Test of matchFormatToPath method, of class RandomAccessGenotypeDataReaderFormats.
	 */
	@Test
	public void testMatchFormatToPath() throws URISyntaxException {
		assertEquals(RandomAccessGenotypeDataReaderFormats.matchFormatToPath(getTriTyperFolder().getAbsolutePath()), RandomAccessGenotypeDataReaderFormats.TRITYPER);
	}

	
}
