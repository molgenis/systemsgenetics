/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.sampleFilter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.ResourceTest;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.trityper.TriTyperGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;
import static org.testng.Assert.*;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

/**
 *
 * @author Patrick Deelen
 */
public class SampleFilterGenotypeDataNGTest extends ResourceTest {
	
	private RandomAccessGenotypeData genotypeDataOriginal;
	private SampleFilterableGenotypeData genotypeDataFiltered;
	
	public SampleFilterGenotypeDataNGTest() {
	}

	@BeforeMethod
	public void setUpMethod() throws Exception {
		
		genotypeDataOriginal = new TriTyperGenotypeData(getTriTyperFolder().getAbsolutePath());
		genotypeDataFiltered = new SampleFilterableGenotypeDataDecorator(genotypeDataOriginal, new SampleIncludedFilter());
		
	}

	/**
	 * Test of getSeqNames method, of class SampleFilterGenotypeData.
	 */
	@Test
	public void testGetSeqNames() {
		assertEquals(genotypeDataFiltered.getSeqNames(), genotypeDataOriginal.getSeqNames());
	}

	/**
	 * Test of getSequences method, of class SampleFilterGenotypeData.
	 */
	@Test
	public void testGetSequences() {
		assertEquals(genotypeDataFiltered.getSequences(), genotypeDataOriginal.getSequences());
	}

	/**
	 * Test of getSequenceByName method, of class SampleFilterGenotypeData.
	 */
	@Test
	public void testGetSequenceByName() {
		assertEquals(genotypeDataFiltered.getSequenceByName("22"), genotypeDataOriginal.getSequenceByName("22"));
	}

	/**
	 * Test of getVariantsByPos method, of class SampleFilterGenotypeData.
	 */
	@Test
	public void testGetVariantsByPos() {
		testFilteredRs11089130(genotypeDataFiltered.getVariantsByPos("22", 14431347).iterator().next());
	}

	/**
	 * Test of getSnpVariantByPos method, of class SampleFilterGenotypeData.
	 */
	@Test
	public void testGetSnpVariantByPos() {
		testFilteredRs11089130(genotypeDataFiltered.getSnpVariantByPos("22", 14431347));
	}

	/**
	 * Test of getSequenceGeneticVariants method, of class SampleFilterGenotypeData.
	 */
	@Test
	public void testGetSequenceGeneticVariants() {
		int counter = 0;
		for(GeneticVariant variant : genotypeDataFiltered.getSequenceGeneticVariants("22")){
			
			if(variant.getPrimaryVariantId().equals("rs11089130")){
				testFilteredRs11089130(variant);				
			}
			
			++counter;
		}
		assertEquals(counter, 9);
	}

	/**
	 * Test of getVariantsByRange method, of class SampleFilterGenotypeData.
	 */
	@Test
	public void testGetVariantsByRange() {
		int counter = 0;
		for(GeneticVariant variant : genotypeDataFiltered.getVariantsByRange("22", 1, 14432619)){
			
			if(variant.getPrimaryVariantId().equals("rs11089130")){
				testFilteredRs11089130(variant);				
			}
			
			++counter;
		}
		assertEquals(counter, 2);
	}

	/**
	 * Test of getVariantAnnotations method, of class SampleFilterGenotypeData.
	 */
	@Test
	public void testGetVariantAnnotations() {
		assertEquals(genotypeDataFiltered.getVariantAnnotations(), genotypeDataOriginal.getVariantAnnotations());
	}

	/**
	 * Test of getVariantAnnotation method, of class SampleFilterGenotypeData.
	 */
	@Test
	public void testGetVariantAnnotation() {
		assertEquals(genotypeDataFiltered.getVariantAnnotations(), genotypeDataOriginal.getVariantAnnotations());
	}

	/**
	 * Test of getSampleAnnotations method, of class SampleFilterGenotypeData.
	 */
	@Test
	public void testGetSampleAnnotations() {
		assertEquals(genotypeDataFiltered.getSampleAnnotations(), genotypeDataOriginal.getSampleAnnotations());
	}

	/**
	 * Test of getSampleAnnotation method, of class SampleFilterGenotypeData.
	 */
	@Test
	public void testGetSampleAnnotation() {
		assertEquals(genotypeDataFiltered.getSampleAnnotation(GenotypeData.BOOL_INCLUDE_SAMPLE), genotypeDataOriginal.getSampleAnnotation(GenotypeData.BOOL_INCLUDE_SAMPLE));
	}

	/**
	 * Test of getSamples method, of class SampleFilterGenotypeData.
	 */
	@Test
	public void testGetSamples() {
		assertEquals(genotypeDataOriginal.getSamples().size(), 9);
		
		boolean encountered_F1047_1047 = false;
		for(Sample sample : genotypeDataOriginal.getSamples()){
			if(sample.getId().equals("F1047-1047")){
				encountered_F1047_1047 = true;
			}
		}
		assertEquals(encountered_F1047_1047, true);
		
		encountered_F1047_1047 = false;
		for(Sample sample : genotypeDataFiltered.getSamples()){
			if(sample.getId().equals("F1047-1047")){
				encountered_F1047_1047 = true;
			}
		}
		assertEquals(encountered_F1047_1047, false);
		
		
		assertEquals(genotypeDataFiltered.getSamples().size(), 8);
	}

	/**
	 * Test of iterator method, of class SampleFilterGenotypeData.
	 */
	@Test
	public void testIterator() {
		int counter = 0;
		for(GeneticVariant variant : genotypeDataFiltered){
			
			if(variant.getPrimaryVariantId().equals("rs11089130")){
				testFilteredRs11089130(variant);				
			}
			
			++counter;
		}
		assertEquals(counter, 10);
	}

	/**
	 * Test of getIncludedSampleCount method, of class SampleFilterGenotypeData.
	 */
	@Test
	public void testGetIncludeCount() {
		assertEquals(genotypeDataFiltered.getIncludedSampleCount(), 8);
	}
	
	private void testFilteredRs11089130(GeneticVariant rs11089130){
		
		List<Alleles> sampleAlleles = rs11089130.getSampleVariants();
		float[] sampleDosage = rs11089130.getSampleDosages();
		byte[] sampleCalledDosage = rs11089130.getSampleCalledDosages();
		List<Boolean> samplePhasing = rs11089130.getSamplePhasing();
		
		assertEquals(sampleAlleles.size(), 8);
		assertEquals(sampleDosage.length, 8);
		assertEquals(sampleCalledDosage.length, 8);
		assertEquals(samplePhasing.size(), 8);
		
		ArrayList<Alleles> expectedAlleles = new ArrayList<Alleles>(8);
		expectedAlleles.add(Alleles.createAlleles(Allele.C, Allele.C));
		expectedAlleles.add(Alleles.createAlleles(Allele.C, Allele.G));
		expectedAlleles.add(Alleles.createAlleles(Allele.G, Allele.G));
		expectedAlleles.add(Alleles.createAlleles(Allele.C, Allele.C));
		expectedAlleles.add(Alleles.createAlleles(Allele.C, Allele.C));
		expectedAlleles.add(Alleles.createAlleles(Allele.C, Allele.C));
		expectedAlleles.add(Alleles.createAlleles(Allele.C, Allele.G));
		expectedAlleles.add(Alleles.createAlleles(Allele.C, Allele.G));
		
		assertEquals(sampleAlleles, expectedAlleles);
		
		float[] expectedDosage = {2, 1, 0, 2, 2, 2, 1, 1};
		
		assertEqualsFloatArray(sampleDosage, expectedDosage);
		
		byte[] expectedCalledDosage = {2, 1, 0, 2, 2, 2, 1, 1};
		
		assertEqualsByteArray(sampleCalledDosage, expectedCalledDosage);
		
				
	}
	
	private void assertEqualsFloatArray(float[] d1, float[] d2){
		
		assertEquals(d1.length, d2.length);
		
		for(int i = 0 ; i < d1.length ; ++i){
			assertEquals(d1[i], d2[i]);
		}
		
	}
	
	private void assertEqualsByteArray(byte[] d1, byte[] d2){
		
		assertEquals(d1.length, d2.length);
		
		for(int i = 0 ; i < d1.length ; ++i){
			assertEquals(d1[i], d2[i]);
		}
		
	}
	
}
