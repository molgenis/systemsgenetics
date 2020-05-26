/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.trityper;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.ResourceTest;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.Sequence;
import org.molgenis.genotype.util.Utils;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variantFilter.VariantQcChecker;
import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertNotNull;
import static org.testng.Assert.assertNull;
import static org.testng.Assert.assertTrue;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

/**
 *
 * @author Patrick Deelen
 */
public class TriTyperGenotypeDataFilteredTest extends ResourceTest {

    
    private TriTyperGenotypeData genotypeData;

	@BeforeClass
	public void beforeClass() throws IOException, URISyntaxException
	{
		genotypeData = new TriTyperGenotypeData(new File(getTriTyperFolder().getAbsolutePath()), 1, new VariantQcChecker(0, 1, 0), true);
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
		System.out.println(variants.toString());
		assertEquals(variants.size(), 6);
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
		assertEquals(sampleVariants.size(), 8);
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
		assertEquals(samples.size(), 8);
		assertEquals(samples.get(0).getId(), "F1042-1042");
		
		boolean foundExcluded = false;
		for(Sample sample : samples){
			if(sample.getId().equals("F1047-1047")){
				foundExcluded = true;
			}
		}
		assertEquals(foundExcluded, false);
		
		assertNull(samples.get(0).getFamilyId());
	}

	@Test
	public void testGetSamplePhasing()
	{
		Iterable <GeneticVariant> variants = genotypeData.getVariantsByPos("22", 14431347);
		
		
		List<Boolean> phasing = genotypeData.getSamplePhasing(variants.iterator().next());
		
		
		assertEquals(phasing,
				Arrays.asList(false, false, false, false, false, false, false, false));
	}

	@Test
	public void testGetSnpVariantByPos()
	{
		int pos = 14433624;
		GeneticVariant variant = genotypeData.getSnpVariantByPos("22", pos);
		assertNotNull(variant);
		assertEquals(variant.getStartPos(), pos);

		ArrayList<Alleles> expectedSampleAlleles = new ArrayList<Alleles>(8);
		expectedSampleAlleles.add(Alleles.createBasedOnChars('G', 'G'));
		expectedSampleAlleles.add(Alleles.createBasedOnChars('A', 'G'));
		expectedSampleAlleles.add(Alleles.createBasedOnChars('G', 'G'));
		expectedSampleAlleles.add(Alleles.createBasedOnChars('G', 'G'));
		expectedSampleAlleles.add(Alleles.createBasedOnChars('G', 'G'));
		expectedSampleAlleles.add(Alleles.createBasedOnChars('G', 'G'));
		expectedSampleAlleles.add(Alleles.createBasedOnChars('G', 'G'));
		expectedSampleAlleles.add(Alleles.createBasedOnChars('A', 'G'));

		assertEquals(variant.getSampleVariants(), expectedSampleAlleles);
		assertEquals(variant.getVariantAlleles(), Alleles.createBasedOnChars('G', 'A'));

		byte[] expectedCalledDosage = new byte[]
		{ 2, 1, 2, 2, 2, 2, 2, 1 };

		assertEquals(variant.getSampleCalledDosages(), expectedCalledDosage);
		
		//Try to get excluded variant
		assertNull(genotypeData.getSnpVariantByPos("22", 14432918));
		
		//now get another variant and do some reloading to test stability. Replaced with other position. 
		GeneticVariant variant2 = genotypeData.getSnpVariantByPos("22", 14431347);
		variant2.getVariantAlleles();
		variant2.getSampleVariants();
		
		assertEquals(variant.getSampleVariants(), expectedSampleAlleles);
		assertEquals(variant.getVariantAlleles(), Alleles.createBasedOnChars('G', 'A'));
		
		GeneticVariant variant4 = genotypeData.getSnpVariantByPos("22", 14433758);
		variant4.getVariantAlleles();
		variant4.getSampleVariants();
		
		GeneticVariant variant3 = genotypeData.getSnpVariantByPos("22", pos);
		assertEquals(variant3.getSampleVariants(), expectedSampleAlleles);
		assertEquals(variant3.getVariantAlleles(), Alleles.createBasedOnChars('G', 'A'));

	}

}
