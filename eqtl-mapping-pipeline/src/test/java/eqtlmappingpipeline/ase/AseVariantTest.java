/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.ase;

import cern.colt.list.tint.IntArrayList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.variant.id.GeneticVariantId;
import static org.testng.Assert.*;
import org.testng.annotations.Test;

/**
 *
 * @author Patrick Deelen
 */
public class AseVariantTest {

	private final AseVariant aseVariant1;
	private final AseVariant aseVariant2;
	private final AseVariant aseVariant3;
	private final AseVariant aseVariant4;

	public AseVariantTest() {
		aseVariant1 = new AseVariant("1", 1, GeneticVariantId.createVariantId("rs1"), Allele.A, Allele.C);
		aseVariant1.addCounts(10, 20, "sample1");
		aseVariant1.addCounts(20, 30, "sample2");
		aseVariant1.addCounts(20, 21, "sample3");
		aseVariant1.addCounts(30, 20, "sample4");

		aseVariant2 = new AseVariant("1", 1, GeneticVariantId.createVariantId("rs1"), Allele.A, Allele.C);
		aseVariant2.addCounts(20, 10, "sample1");

		aseVariant3 = new AseVariant("1", 1, GeneticVariantId.createVariantId("rs1"), Allele.A, Allele.C);
		aseVariant3.addCounts(124, 99, "sample1");
		aseVariant3.addCounts(28, 179, "sample2");

		aseVariant4 = new AseVariant("1", 1, GeneticVariantId.createVariantId("rs1"), Allele.A, Allele.C);
		aseVariant4.addCounts(211, 27, "sample1");
		aseVariant4.addCounts(196, 45, "sample2");
		aseVariant4.addCounts(187, 54, "sample3");
		
	}

	/**
	 * Test of getChr method, of class AseVariant.
	 */
	@Test
	public void testGetChr() {
		assertEquals(aseVariant1.getChr(), "1");
	}

	/**
	 * Test of getPos method, of class AseVariant.
	 */
	@Test
	public void testGetPos() {
		assertEquals(aseVariant1.getPos(), 1);
	}

	/**
	 * Test of getId method, of class AseVariant.
	 */
	@Test
	public void testGetId() {
		assertEquals(aseVariant1.getId().getPrimairyId(), "rs1");
	}

	/**
	 * Test of getA1 method, of class AseVariant.
	 */
	@Test
	public void testGetA1() {
		assertEquals(aseVariant1.getA1(), Allele.A);
	}

	/**
	 * Test of getA2 method, of class AseVariant.
	 */
	@Test
	public void testGetA2() {
		assertEquals(aseVariant1.getA2(), Allele.C);
	}

	/**
	 * Test of getA1Counts method, of class AseVariant.
	 */
	@Test
	public void testGetA1Counts() {
		IntArrayList expResult = new IntArrayList(new int[]{10, 20, 20, 30});
		assertEquals(aseVariant1.getA1Counts(), expResult);
	}

	/**
	 * Test of getA2Counts method, of class AseVariant.
	 */
	@Test
	public void testGetA2Counts() {
		IntArrayList expResult = new IntArrayList(new int[]{20, 30, 21, 20});
		assertEquals(aseVariant1.getA2Counts(), expResult);
	}

	@Test
	public void testGetZscore(){
		assertEquals(aseVariant1.getMetaZscore(), 0.8255035, 0.0001);
		assertEquals(aseVariant2.getMetaZscore(), -1.651007, 0.0001);
		assertEquals(aseVariant3.getMetaZscore(), 6.639489, 0.00001);
		
	}

	@Test
	public void testGetPvalues(){
		assertEquals(aseVariant1.getMetaPvalue(), 0.4090858, 0.00001);
		assertEquals(aseVariant2.getMetaPvalue(), 0.09873715, 0.00001);
		assertEquals(aseVariant3.getMetaPvalue(), 3.147722e-11, 0.00001);
	}

	@Test
	public void testGetCountPearsonR(){

		assertEquals(aseVariant4.getCountPearsonR(), -0.9989061, 0.00001);

	}
	
	@Test
	public void testGetSampleIds(){
		List<String> samples = Arrays.asList(new String[]{"sample1", "sample2", "sample3", "sample4"});
		assertEquals(aseVariant1.getSampleIds(), samples);
		
	}

}