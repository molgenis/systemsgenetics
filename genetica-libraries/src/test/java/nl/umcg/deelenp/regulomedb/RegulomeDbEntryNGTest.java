/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.deelenp.regulomedb;

import umcg.genetica.io.regulomedb.RegulomeDbEntry;
import umcg.genetica.io.regulomedb.RegulomeDbSupportingData;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import static org.testng.Assert.*;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

/**
 *
 * @author Patrick Deelen
 */
public class RegulomeDbEntryNGTest {
	
	public RegulomeDbEntryNGTest() {
	}
	
	private RegulomeDbEntry instance;

	@BeforeMethod
	public void setUpMethod() throws Exception {
		instance = new RegulomeDbEntry("chr1\t20206\trs12239663\tMotifs|PWM|Zic1, Motifs|PWM|Zic3, Chromatin_Structure|DNase-seq|Sknmc\t5");
	}

	/**
	 * Test of getChr method, of class RegulomeDbEntry.
	 */
	@Test
	public void testGetChr() {
		System.out.println("getChr");
		String expResult = "chr1";
		String result = instance.getChr();
		assertEquals(result, expResult);
	}

	/**
	 * Test of getChrPos method, of class RegulomeDbEntry.
	 */
	@Test
	public void testGetChrPos() {
		System.out.println("getChrPos");
		int expResult = 20206;
		int result = instance.getChrPos();
		assertEquals(result, expResult);
	}

	/**
	 * Test of getVariantId method, of class RegulomeDbEntry.
	 */
	@Test
	public void testGetVariantId() {
		System.out.println("getVariantId");
		String expResult = "rs12239663";
		String result = instance.getVariantId();
		assertEquals(result, expResult);
	}

	/**
	 * Test of getRegulomeDbScore method, of class RegulomeDbEntry.
	 */
	@Test
	public void testGetRegulomeDbScore() {
		System.out.println("getRegulomeDbScore");
		String expResult = "5";
		String result = instance.getRegulomeDbScore();
		assertEquals(result, expResult);
	}

	/**
	 * Test of getSupportData method, of class RegulomeDbEntry.
	 */
	@Test
	public void testGetSupportData() {
		System.out.println("getSupportData");
		Map<String, List<RegulomeDbSupportingData>> instanceSupportData = instance.getSupportData();
		
		assertEquals(instanceSupportData.containsKey("Motifs"), true);
		assertEquals(instanceSupportData.containsKey("Chromatin_Structure"), true);
		assertEquals(instanceSupportData.containsKey("Protein_Binding"), false);
		
		assertEquals(instanceSupportData.get("Motifs").size(), 2);
		assertEquals(instanceSupportData.get("Chromatin_Structure").size(), 1);
		
		Iterator<RegulomeDbSupportingData> suppporintDataIterator = instanceSupportData.get("Motifs").iterator();
		RegulomeDbSupportingData supportingDataEntry = suppporintDataIterator.next();
		assertEquals(supportingDataEntry.getSupportClass(), "Motifs");
		assertEquals(supportingDataEntry.getSupportMethod(), "PWM");
		assertEquals(supportingDataEntry.getSupportValue(), "Zic1");

		supportingDataEntry = suppporintDataIterator.next();
		assertEquals(supportingDataEntry.getSupportClass(), "Motifs");
		assertEquals(supportingDataEntry.getSupportMethod(), "PWM");
		assertEquals(supportingDataEntry.getSupportValue(), "Zic3");
		
		suppporintDataIterator = instanceSupportData.get("Chromatin_Structure").iterator();
		supportingDataEntry = suppporintDataIterator.next();
		assertEquals(supportingDataEntry.getSupportClass(), "Chromatin_Structure");
		assertEquals(supportingDataEntry.getSupportMethod(), "DNase-seq");
		assertEquals(supportingDataEntry.getSupportValue(), "Sknmc");
		
	}

	/**
	 * Test of hashCode method, of class RegulomeDbEntry.
	 */
	@Test
	public void testHashCode() throws Exception {
		System.out.println("hashCode");
		RegulomeDbEntry other = new RegulomeDbEntry("chr1\t20206\trs12239663\tMotifs|PWM|Zic1, Motifs|PWM|Zic3, Chromatin_Structure|DNase-seq|Sknmc\t5");
		int expResult = other.hashCode();
		int result = instance.hashCode();
		assertEquals(result, expResult);
	}

	/**
	 * Test of equals method, of class RegulomeDbEntry.
	 */
	@Test
	public void testEquals() throws Exception {
		System.out.println("equals");
		RegulomeDbEntry other = new RegulomeDbEntry("chr1\t20206\trs12239663\tMotifs|PWM|Zic1, Motifs|PWM|Zic3, Chromatin_Structure|DNase-seq|Sknmc\t5");
		boolean expResult = true;
		boolean result = instance.equals(other);
		assertEquals(result, expResult);
		
		
		other = new RegulomeDbEntry("chr2\t20206\trs12239663\tMotifs|PWM|Zic1, Motifs|PWM|Zic3, Chromatin_Structure|DNase-seq|Sknmc\t5");
		expResult = false;
		result = instance.equals(other);
		assertEquals(result, expResult);
		
		other = new RegulomeDbEntry("chr1\t10206\trs12239663\tMotifs|PWM|Zic1, Motifs|PWM|Zic3, Chromatin_Structure|DNase-seq|Sknmc\t5");
		expResult = false;
		result = instance.equals(other);
		assertEquals(result, expResult);
		
		other = new RegulomeDbEntry("chr1\t20206\t12239663\tMotifs|PWM|Zic1, Motifs|PWM|Zic3, Chromatin_Structure|DNase-seq|Sknmc\t5");
		expResult = false;
		result = instance.equals(other);
		assertEquals(result, expResult);
		
		other = new RegulomeDbEntry("chr1\t20206\trs12239663\tM|PWM|Zic1, Motifs|PWM|Zic3, Chromatin_Structure|DNase-seq|Sknmc\t5");
		expResult = false;
		result = instance.equals(other);
		assertEquals(result, expResult);
		
		other = new RegulomeDbEntry("chr1\t20206\trs12239663\tMotifs|PWM|Zic1, Motifs|PWM|Zic3, Chromatin_Structure|DNase-seq|Sknmc\t6");
		expResult = false;
		result = instance.equals(other);
		assertEquals(result, expResult);
		
		other = new RegulomeDbEntry("chr1\t20206\trs12239663\tMotifs|PWM|Zic1, Motifs|PWM|Zic3, Chromatin_Structure|DNase-seq|Bazinga\t5");
		expResult = false;
		result = instance.equals(other);
		assertEquals(result, expResult);
		
		other = new RegulomeDbEntry("chr1\t20206\trs12239663\tMotifs|PWM|Zic1, Motifs|PWM|Zic, Chromatin_Structure|DNase-seq|Sknmc\t5");
		expResult = false;
		result = instance.equals(other);
		assertEquals(result, expResult);
		
	}

	/**
	 * Test of hasSupportDataClass method, of class RegulomeDbEntry.
	 */
	@Test
	public void testHasSupportDataClass() {
		System.out.println("hasSupportDataClass");
		
		assertEquals(instance.hasSupportDataClass("Motifs"), true);
		assertEquals(instance.hasSupportDataClass("Chromatin_Structure"), true);
		assertEquals(instance.hasSupportDataClass("Protein_Binding"), false);

	}
}
