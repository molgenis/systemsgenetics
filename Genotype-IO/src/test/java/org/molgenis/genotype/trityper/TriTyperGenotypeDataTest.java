/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.trityper;

import java.io.IOException;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.ResourceTest;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.Sequence;
import org.molgenis.genotype.util.Utils;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variantFilter.VariantFilterBiAllelic;
import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertNotNull;
import static org.testng.Assert.assertNull;
import static org.testng.Assert.assertTrue;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

/**
 *
 * @author harmjan
 */
public class TriTyperGenotypeDataTest extends ResourceTest {
    
    
    private TriTyperGenotypeData genotypeData;

	@BeforeClass
	public void beforeClass() throws IOException, URISyntaxException
	{
		genotypeData = new TriTyperGenotypeData(getTriTyperFolder().getAbsolutePath());
	}
    
    @Test
	public void getSeqNames()
	{
		List<String> seqNames = genotypeData.getSeqNames();
		assertNotNull(seqNames);
		assertEquals(seqNames.size(), 2);
		assertEquals(seqNames.get(0), "22");
		assertEquals(seqNames.get(1), "23");
	}

	@Test
	public void testGetSequences()
	{
		Iterable<Sequence> sequences = genotypeData.getSequences();
		assertNotNull(sequences);
        
        int count = 0;
        for(Sequence s : sequences){
            ++count;
        }
        
		assertEquals(count, 2);
	}

	@Test
	public void testGetSequenceByName() throws IOException
	{
		Sequence sequence = genotypeData.getSequenceByName("22");
		assertNotNull(sequence);
		assertEquals(sequence.getName(), "22");

		List<GeneticVariant> variants = Utils.iteratorToList(genotypeData.getSequenceGeneticVariants(sequence.getName()).iterator());
		assertNotNull(variants);
		assertEquals(variants.size(), 9);
		GeneticVariant variant = variants.get(0);
		assertEquals(variant.getPrimaryVariantId(), "rs11089130");
		assertEquals(variant.getStartPos(), 14431347);

		assertEquals(variant.getSequenceName(), "22");
		assertEquals(variant.isSnp(), true);

		List<String> alleles = variant.getVariantAlleles().getAllelesAsString();
		assertNotNull(alleles);
		assertEquals(alleles.size(), 2);
		assertTrue(alleles.contains("C"));
		assertTrue(alleles.contains("G"));

		List<Alleles> sampleVariants = variant.getSampleVariants();
		assertNotNull(sampleVariants);
		assertEquals(sampleVariants.size(), 9);
		assertNotNull(sampleVariants.get(0).getAllelesAsChars());
		assertEquals(sampleVariants.get(0).getAlleles().size(), 2);
		assertEquals(sampleVariants.get(0).getAllelesAsChars()[0], 'C');
		assertEquals(sampleVariants.get(0).getAllelesAsChars()[1], 'C');
		assertNotNull(sampleVariants.get(1).getAllelesAsChars());
		assertEquals(sampleVariants.get(1).getAllelesAsChars().length, 2);
		assertEquals(sampleVariants.get(1).getAllelesAsChars()[0], 'C');
		assertEquals(sampleVariants.get(1).getAllelesAsChars()[1], 'G');
		assertNotNull(sampleVariants.get(2).getAllelesAsChars());
		assertEquals(sampleVariants.get(2).getAllelesAsChars().length, 2);
		assertEquals(sampleVariants.get(2).getAllelesAsChars()[0], 'G');
		assertEquals(sampleVariants.get(2).getAllelesAsChars()[1], 'G');
	}

	@Test
	public void testGetSamples()
	{
		List<Sample> samples = genotypeData.getSamples();
		assertNotNull(samples);
		assertEquals(samples.size(), 9);
		assertEquals(samples.get(0).getId(), "F1042-1042");
		assertNull(samples.get(0).getFamilyId());
	}

	@Test
	public void testGetSamplePhasing()
	{
		Iterable <GeneticVariant> variants = genotypeData.getVariantsByPos("22", 14431347);
		
		
		List<Boolean> phasing = genotypeData.getSamplePhasing(variants.iterator().next());
		
		
		assertEquals(phasing,
				Arrays.asList(false, false, false, false, false, false, false, false, false));
	}
	
	/**
	 * Test of iterator method, of class SampleFilterGenotypeData.
	 */
	@Test
	public void testIterator() {
		int counter = 0;
		for(GeneticVariant variant : genotypeData){
			
			if(variant.getPrimaryVariantId().equals("rs11089130")){
				testFilteredRs11089130(variant);				
			}
			
			++counter;
		}
		assertEquals(counter, 10);
	}

	@Test
	public void testGetSnpVariantByPos()
	{
		int pos = 14433624;
		GeneticVariant variant = genotypeData.getSnpVariantByPos("22", pos);
		assertNotNull(variant);
		assertEquals(variant.getStartPos(), pos);

		ArrayList<Alleles> expectedSampleAlleles = new ArrayList<Alleles>(9);
		expectedSampleAlleles.add(Alleles.createBasedOnChars('G', 'G'));
		expectedSampleAlleles.add(Alleles.createBasedOnChars('A', 'G'));
		expectedSampleAlleles.add(Alleles.createBasedOnChars('G', 'G'));
		expectedSampleAlleles.add(Alleles.createBasedOnChars('G', 'G'));
		expectedSampleAlleles.add(Alleles.createBasedOnChars('G', 'G'));
		expectedSampleAlleles.add(Alleles.createBasedOnChars('A', 'G'));
		expectedSampleAlleles.add(Alleles.createBasedOnChars('G', 'G'));
		expectedSampleAlleles.add(Alleles.createBasedOnChars('G', 'G'));
		expectedSampleAlleles.add(Alleles.createBasedOnChars('A', 'G'));

		assertEquals(variant.getSampleVariants(), expectedSampleAlleles);
		assertEquals(variant.getVariantAlleles(), Alleles.createBasedOnChars('G', 'A'));

		byte[] expectedCalledDosage = new byte[]
		{ 2, 1, 2, 2, 2, 1, 2, 2, 1 };

		assertEquals(variant.getSampleCalledDosages(), expectedCalledDosage);

	}
	
	@Test
	public void getVariantByIdTest(){
		
		GeneticVariant variant = genotypeData.getVariantIdMap().get("rs7510853");
		
		assertEquals(variant.getPrimaryVariantId(), "rs7510853");
		
		assertEquals(variant.getVariantAlleles(), Alleles.createAlleles(Allele.C));
		
		
	}

	@Test
	public void filteredVariantIdMapTest(){
		assertNull(genotypeData.getVariantIdMap(new VariantFilterBiAllelic()).get("rs7510853"));
	}
	
	private void testFilteredRs11089130(GeneticVariant rs11089130){
		
		List<Alleles> sampleAlleles = rs11089130.getSampleVariants();
		float[] sampleDosage = rs11089130.getSampleDosages();
		byte[] sampleCalledDosage = rs11089130.getSampleCalledDosages();
		List<Boolean> samplePhasing = rs11089130.getSamplePhasing();
		
		assertEquals(sampleAlleles.size(), 9);
		assertEquals(sampleDosage.length, 9);
		assertEquals(sampleCalledDosage.length, 9);
		assertEquals(samplePhasing.size(), 9);
		
		ArrayList<Alleles> expectedAlleles = new ArrayList<Alleles>(9);
		expectedAlleles.add(Alleles.createAlleles(Allele.C, Allele.C));
		expectedAlleles.add(Alleles.createAlleles(Allele.C, Allele.G));
		expectedAlleles.add(Alleles.createAlleles(Allele.G, Allele.G));
		expectedAlleles.add(Alleles.createAlleles(Allele.C, Allele.C));
		expectedAlleles.add(Alleles.createAlleles(Allele.C, Allele.C));
		expectedAlleles.add(Alleles.createAlleles(Allele.C, Allele.C));
		expectedAlleles.add(Alleles.createAlleles(Allele.C, Allele.C));
		expectedAlleles.add(Alleles.createAlleles(Allele.C, Allele.G));
		expectedAlleles.add(Alleles.createAlleles(Allele.C, Allele.G));
		
		assertEquals(sampleAlleles, expectedAlleles);
		
		float[] expectedDosage = {2, 1, 0, 2, 2, 2, 2, 1, 1};
		
		assertEqualsFloatArray(sampleDosage, expectedDosage);
		
		byte[] expectedCalledDosage = {2, 1, 0, 2, 2, 2, 2, 1, 1};
		
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
