/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.ase;

import cern.colt.list.tdouble.DoubleArrayList;
import cern.colt.list.tint.IntArrayList;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.variant.id.GeneticVariantId;
import static org.testng.Assert.*;
import org.testng.annotations.Test;

/**
 *
 * @author Patrick Deelen
 */
public class AseVariantTest {
	
	AseVariant aseVariant;
	
	public AseVariantTest() {
		aseVariant = new AseVariant("1", 1, GeneticVariantId.createVariantId("rs1"), Allele.A, Allele.C);
		aseVariant.addCounts(10, 20);
		aseVariant.addCounts(20, 30);
		aseVariant.addCounts(20, 21);
		aseVariant.addCounts(30, 20);
	}

	/**
	 * Test of getChr method, of class AseVariant.
	 */
	@Test
	public void testGetChr() {
		assertEquals(aseVariant.getChr(), "1");
	}

	/**
	 * Test of getPos method, of class AseVariant.
	 */
	@Test
	public void testGetPos() {
		assertEquals(aseVariant.getPos(), 1);
	}

	/**
	 * Test of getId method, of class AseVariant.
	 */
	@Test
	public void testGetId() {
		assertEquals(aseVariant.getId().getPrimairyId(), "rs1");
	}

	/**
	 * Test of getA1 method, of class AseVariant.
	 */
	@Test
	public void testGetA1() {
		assertEquals(aseVariant.getA1(), Allele.A);
	}

	/**
	 * Test of getA2 method, of class AseVariant.
	 */
	@Test
	public void testGetA2() {
		assertEquals(aseVariant.getA2(), Allele.C);
	}

	/**
	 * Test of getA1Counts method, of class AseVariant.
	 */
	@Test
	public void testGetA1Counts() {
		IntArrayList expResult = new IntArrayList(new int[]{10, 20, 20, 30});
		assertEquals(aseVariant.getA1Counts(), expResult);
	}

	/**
	 * Test of getA2Counts method, of class AseVariant.
	 */
	@Test
	public void testGetA2Counts() {
		IntArrayList expResult = new IntArrayList(new int[]{20, 30, 21, 20});
		assertEquals(aseVariant.getA2Counts(), expResult);
	}
	
	public void testGetZscore(){
		assertEquals(aseVariant.getMetaZscore(), 0.8255035, 0.0001);
	}

}