/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.deelenp.regulomedb;

import umcg.genetica.io.regulomedb.RegulomeDbSupportingData;
import static org.testng.Assert.*;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

/**
 *
 * @author Patrick Deelen
 */
public class RegulomeDbSupportingDataNGTest {
	
	private static RegulomeDbSupportingData testData;
	
	public RegulomeDbSupportingDataNGTest() {
	}

	@BeforeClass
	public static void setUpClass() throws Exception {
		testData = new RegulomeDbSupportingData("Chromatin_Structure|DNase-seq|Nb4");
	}

	/**
	 * Test of getSupportClass method, of class RegulomeDbSupportingData.
	 */
	@Test
	public void testGetSupportClass() {
		assertEquals(testData.getSupportClass(), "Chromatin_Structure");
	}

	/**
	 * Test of getSupportMethod method, of class RegulomeDbSupportingData.
	 */
	@Test
	public void testGetSupportMethod() {
		assertEquals(testData.getSupportMethod(), "DNase-seq");
	}

	/**
	 * Test of getSupportValue method, of class RegulomeDbSupportingData.
	 */
	@Test
	public void testGetSupportValue() {
		assertEquals(testData.getSupportValue(), "Nb4");
	}
	
	@Test
	public void testSpecialCase1() throws Exception{
		RegulomeDbSupportingData tmp = new RegulomeDbSupportingData("Chromatin_Structure|DNase-seq");
		assertEquals(tmp.getSupportClass(), "Chromatin_Structure");
		assertEquals(tmp.getSupportMethod(), "DNase-seq");
		assertEquals(tmp.getSupportValue(), "");
	}
	
	@Test
	public void testSpecialCase2() throws Exception{
		RegulomeDbSupportingData tmp = new RegulomeDbSupportingData("Chromatin_Structure|DNase-seq|");
		assertEquals(tmp.getSupportClass(), "Chromatin_Structure");
		assertEquals(tmp.getSupportMethod(), "DNase-seq");
		assertEquals(tmp.getSupportValue(), "");
	}
	
	@Test
	public void testSpecialCase3() throws Exception{
		RegulomeDbSupportingData tmp = new RegulomeDbSupportingData("Chromatin_Structure|DNase-seq|a|b");
		assertEquals(tmp.getSupportClass(), "Chromatin_Structure");
		assertEquals(tmp.getSupportMethod(), "DNase-seq");
		assertEquals(tmp.getSupportValue(), "a|b");
	}

	@Test(expectedExceptions = Exception.class)
	public void testExceptionTooLongSupportingData() throws Exception{
		RegulomeDbSupportingData tmp = new RegulomeDbSupportingData("Chromatin_Structure|DNase-seq|Nb4|Bazinga|test");
	}

	/**
	 * Test of hashCode method, of class RegulomeDbSupportingData.
	 */
	@Test
	public void testHashCode() throws Exception {
		System.out.println("hashCode");
		RegulomeDbSupportingData instance = testData;
		RegulomeDbSupportingData otherData = new RegulomeDbSupportingData("Chromatin_Structure|DNase-seq|Nb4");
		int expResult = otherData.hashCode();
		int result = instance.hashCode();
		assertEquals(result, expResult);
	}

	/**
	 * Test of equals method, of class RegulomeDbSupportingData.
	 */
	@Test
	public void testEquals() throws Exception {
		System.out.println("equals");
		RegulomeDbSupportingData instance = testData;
		RegulomeDbSupportingData otherData = new RegulomeDbSupportingData("Chromatin_Structure|DNase-seq|Nb4");
		boolean expResult = true;
		boolean result = instance.equals(otherData);
		assertEquals(result, expResult);
		
		otherData = new RegulomeDbSupportingData("Chromatin_Structure|DNase-seq|Nb3");
		expResult = false;
		result = instance.equals(otherData);
		assertEquals(result, expResult);
		
		otherData = new RegulomeDbSupportingData("Chromatin_Structure|DNase-seq-Bazinga|Nb4");
		expResult = false;
		result = instance.equals(otherData);
		assertEquals(result, expResult);
		
		otherData = new RegulomeDbSupportingData("CS|DNase-seq|Nb4");
		expResult = false;
		result = instance.equals(otherData);
		assertEquals(result, expResult);
		
	}

	
}
