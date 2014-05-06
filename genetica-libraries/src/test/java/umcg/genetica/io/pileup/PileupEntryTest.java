/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.pileup;

import umcg.genetica.io.pileup.PileupEntry;
import umcg.genetica.io.pileup.PileupParseException;
import org.molgenis.genotype.Allele;
import static org.testng.Assert.*;
import org.testng.annotations.Test;

/**
 *
 * @author Patrick Deelen
 */
public class PileupEntryTest {
	
	PileupEntry entry1;
	PileupEntry entry2;
	PileupEntry entry3;
	PileupEntry entry4;
	
	public PileupEntryTest() throws PileupParseException {
		
		entry1 = new PileupEntry("1", 1, Allele.G, 10, "...........A");
		entry2 = new PileupEntry("1", 1, Allele.G, 10, ".$...........^~.");
		entry3 = new PileupEntry("1", 1, Allele.A, 10, ".$......+2AG.+2AG.+2AGGG");
		entry4 = new PileupEntry("1", 1, Allele.A, 10, ".,AGCTNagctn^~$<>-3ACG+12AAACCCTTTGGG");
		
	}

	/**
	 * Test of getChr method, of class PileupEntry.
	 */
	@Test
	public void testGetChr() {
		assertEquals(entry1.getChr(), "1");
	}

	/**
	 * Test of getPos method, of class PileupEntry.
	 */
	@Test
	public void testGetPos() {
		assertEquals(entry1.getPos(), 1);
	}

	/**
	 * Test of getRefAllele method, of class PileupEntry.
	 */
	@Test
	public void testGetRefAllele() {
		assertEquals(entry1.getRefAllele(), Allele.G);
	}

	/**
	 * Test of getReadDepth method, of class PileupEntry.
	 */
	@Test
	public void testGetReadDepth() {
		assertEquals(entry1.getReadDepth(), 10);
	}

	/**
	 * Test of getAlleleCounts method, of class PileupEntry.
	 */
	@Test
	public void testGetAlleleCounts() {
		
		assertEquals(entry1.getAlleleCounts().get(Allele.A), 1);
		assertEquals(entry1.getAlleleCounts().get(Allele.C), 0);
		assertEquals(entry1.getAlleleCounts().get(Allele.G), 11);
		assertEquals(entry1.getAlleleCounts().get(Allele.T), 0);
		
		assertEquals(entry2.getAlleleCounts().get(Allele.A), 0);
		assertEquals(entry2.getAlleleCounts().get(Allele.C), 0);
		assertEquals(entry2.getAlleleCounts().get(Allele.G), 13);
		assertEquals(entry2.getAlleleCounts().get(Allele.T), 0);

		assertEquals(entry3.getAlleleCounts().get(Allele.A), 9);
		assertEquals(entry3.getAlleleCounts().get(Allele.C), 0);
		assertEquals(entry3.getAlleleCounts().get(Allele.G), 2);
		assertEquals(entry3.getAlleleCounts().get(Allele.T), 0);
		
		assertEquals(entry4.getAlleleCounts().get(Allele.A), 4);
		assertEquals(entry4.getAlleleCounts().get(Allele.C), 2);
		assertEquals(entry4.getAlleleCounts().get(Allele.G), 2);
		assertEquals(entry4.getAlleleCounts().get(Allele.T), 2);
	
	}
}