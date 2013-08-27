/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.variantFilter;

import java.util.HashSet;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.ResourceTest;
import org.molgenis.genotype.variant.GeneticVariant;
import static org.testng.Assert.*;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

/**
 *
 * @author Patrick Deelen
 */
public class VariantFilterableGenotypeDataDecoratorNGTest extends ResourceTest {

	private RandomAccessGenotypeData genotypeData;
	HashSet<String> snpsToInclude;
	
	public VariantFilterableGenotypeDataDecoratorNGTest() {
	}

	@BeforeMethod
	public void setUpMethod() throws Exception {

		String testPedPath = getTestPed().toString();
		String pedMapBasePath = testPedPath.substring(0, testPedPath.length() - 4);
		
		snpsToInclude = new HashSet<String>();
		snpsToInclude.add("rs7510853");
		snpsToInclude.add("rs9604721");

		VariantFilter filter = new VariantIdIncludeFilter(snpsToInclude);
		
		genotypeData = RandomAccessGenotypeDataReaderFormats.PED_MAP.createFilteredGenotypeData(pedMapBasePath, 100, filter, null);
		
		
	}

	/**
	 * Test of getSeqNames method, of class
	 * VariantFilterableGenotypeDataDecorator.
	 */
	@Test
	public void testGetSeqNames() {
		assertEquals(genotypeData.getSeqNames().size(), 2);
		assertEquals(genotypeData.getSeqNames().get(0), "22");
	}

	/**
	 * Test of getSequences method, of class
	 * VariantFilterableGenotypeDataDecorator.
	 */
	@Test
	public void testGetSequences() {
		assertEquals(genotypeData.getSequences().iterator().next().getName(), "22");
	}

	/**
	 * Test of getSequenceByName method, of class
	 * VariantFilterableGenotypeDataDecorator.
	 */
	@Test
	public void testGetSequenceByName() {
		assertNotNull(genotypeData.getSequenceByName("22"));
		assertNull(genotypeData.getSequenceByName("1"));
	}

	/**
	 * Test of getVariantsByPos method, of class
	 * VariantFilterableGenotypeDataDecorator.
	 */
	@Test
	public void testGetVariantsByPos() {
		assertEquals(genotypeData.getVariantsByPos("22", 14433624).iterator().hasNext(), false);
		
		assertEquals(genotypeData.getVariantsByPos("22", 14432918).iterator().next().getPrimaryVariantId(), "rs7510853");
		
	}

	/**
	 * Test of getSnpVariantByPos method, of class
	 * VariantFilterableGenotypeDataDecorator.
	 */
	@Test
	public void testGetSnpVariantByPos() {
		assertNull(genotypeData.getSnpVariantByPos("22", 14433624));
		assertEquals(genotypeData.getSnpVariantByPos("22", 14432918).getPrimaryVariantId(), "rs7510853");
	}

	/**
	 * Test of getSequenceGeneticVariants method, of class
	 * VariantFilterableGenotypeDataDecorator.
	 */
	@Test
	public void testGetSequenceGeneticVariants() {
		int count = 0;
		
		for(GeneticVariant variant : genotypeData.getSequenceGeneticVariants("22")){
			assertEquals(snpsToInclude.contains(variant.getPrimaryVariantId()), true);
			++count;
		}
		
		assertEquals(count, 2);
	}

	/**
	 * Test of getVariantsByRange method, of class
	 * VariantFilterableGenotypeDataDecorator.
	 */
	@Test
	public void testGetVariantsByRange() {
		
		int count = 0;
		
		for(GeneticVariant variant : genotypeData.getVariantsByRange("23", 0, Integer.MAX_VALUE)){
			assertEquals(snpsToInclude.contains(variant.getPrimaryVariantId()), true);
			++count;
		}
		
		assertEquals(count, 0);
		
		count = 0;
		
		for(GeneticVariant variant : genotypeData.getVariantsByRange("22", 0, Integer.MAX_VALUE)){
			assertEquals(snpsToInclude.contains(variant.getPrimaryVariantId()), true);
			++count;
		}
		
		assertEquals(count, snpsToInclude.size());
		
	}


	/**
	 * Test of getSamples method, of class
	 * VariantFilterableGenotypeDataDecorator.
	 */
	@Test
	public void testGetSamples() {
		assertEquals(genotypeData.getSamples().size(), 9);
	}

	/**
	 * Test of iterator method, of class VariantFilterableGenotypeDataDecorator.
	 */
	@Test
	public void testIterator() {
		
		int count = 0;
		
		for(GeneticVariant variant : genotypeData){
			assertEquals(snpsToInclude.contains(variant.getPrimaryVariantId()), true);
			++count;
		}
		
		assertEquals(count, snpsToInclude.size());
		
	}

}
