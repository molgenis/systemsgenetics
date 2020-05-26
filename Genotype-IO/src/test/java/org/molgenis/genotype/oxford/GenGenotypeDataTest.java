/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.oxford;

import java.io.IOException;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.ResourceTest;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.Sequence;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.util.Utils;
import org.molgenis.genotype.variant.GeneticVariant;
import static org.molgenis.genotype.util.AssertExtended.*;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

/**
 *
 * @author Patrick Deelen
 */
public class GenGenotypeDataTest extends ResourceTest {

	private GenGenotypeData genotypeData;

	@BeforeClass
	public void beforeClass() throws IOException, URISyntaxException {
		genotypeData = new GenGenotypeData(getTest2Gen(), getTest2Sample());
	}

	public GenGenotypeDataTest() {
	}

	/**
	 * Test of getSequences method, of class GenGenotypeData.
	 */
	@Test
	public void testGetSequences() {
		List<Sequence> sequences = genotypeData.getSequences();
		assertNotNull(sequences);
		assertEquals(sequences.size(), 3);
	}

	@Test
	public void testGetSequenceByName() throws IOException {
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

	/**
	 * Test of getSamples method, of class GenGenotypeData.
	 */
	@Test
	public void testGetSamples() {
		List<Sample> samples = genotypeData.getSamples();
		assertNotNull(samples);
		assertEquals(samples.size(), 9);
		assertEquals(samples.get(0).getId(), "1042");
		assertEquals(samples.get(0).getFamilyId(), "F1042");
		
		assertEquals(samples.get(8).getMissingRate(), 0.1f);
		assertEquals(samples.get(0).getMissingRate(), 0f);
		
		assertEquals((Float) samples.get(0).getAnnotationValues().get("phenotype"), 0f);
		assertEquals((Float) samples.get(1).getAnnotationValues().get("phenotype"), 1f);
		assertEquals((Float) samples.get(2).getAnnotationValues().get("phenotype"), 2f);
		assertTrue(Float.isNaN((Float) samples.get(3).getAnnotationValues().get("phenotype")));
		
	}

	/**
	 * Test of getSampleAnnotationsMap method, of class GenGenotypeData.
	 */
	@Test
	public void testGetSampleAnnotationsMap() {
		
		Map<String, SampleAnnotation> annotationMap = genotypeData.getSampleAnnotationsMap();
		
		int i = 0;
		for(Map.Entry<String, SampleAnnotation> annotation : annotationMap.entrySet()){
			switch(i){
				case 0:
					assertEquals(annotation.getKey(), GenotypeData.SAMPLE_MISSING_RATE_FLOAT);
					assertEquals(annotation.getValue().getName(), "missing");
					assertEquals(annotation.getValue().getSampleAnnotationType(), SampleAnnotation.SampleAnnotationType.OTHER);
					assertEquals(annotation.getValue().getType(), Annotation.Type.FLOAT);
					assertFalse(annotation.getValue().isList());
					break;
				case 1:
					assertEquals(annotation.getKey(), "sex");
					assertEquals(annotation.getValue().getName(), "sex");
					assertEquals(annotation.getValue().getSampleAnnotationType(), SampleAnnotation.SampleAnnotationType.COVARIATE);
					assertEquals(annotation.getValue().getType(), Annotation.Type.STRING);
					assertFalse(annotation.getValue().isList());
					break;
					
				case 2:
					assertEquals(annotation.getKey(), "phenotype");
					assertEquals(annotation.getValue().getName(), "phenotype");
					assertEquals(annotation.getValue().getSampleAnnotationType(), SampleAnnotation.SampleAnnotationType.PHENOTYPE);
					assertEquals(annotation.getValue().getType(), Annotation.Type.FLOAT);
					assertFalse(annotation.getValue().isList());
					break;
					
				default:
					fail("Error expected only 3 annotations");
						
			}
			
			++i;
		}
		
		
	}

	/**
	 * Test of getSamplePhasing method, of class GenGenotypeData.
	 */
	@Test
	public void testGetSamplePhasing() {
		Iterable<GeneticVariant> variants = genotypeData.getVariantsByPos("22", 14431347);

		int counted = 0;
		for (GeneticVariant variant : variants) {

			assertEquals(variant.getSamplePhasing(), Arrays.asList(false, false, false, false, false, false, false, false, false));
			++counted;
		}

		assertEquals(counted, 1);
	}

	/**
	 * Test of getSeqNames method, of class GenGenotypeData.
	 */
	@Test
	public void testGetSeqNames() {
		List<String> seqNames = genotypeData.getSeqNames();
		assertNotNull(seqNames);
		assertEquals(seqNames.size(), 3);
		assertEquals(seqNames.get(0), "22");
		assertEquals(seqNames.get(1), "23");
                assertEquals(seqNames.get(2), "1");
	}

	/**
	 * Test of isOnlyContaingSaveProbabilityGenotypes method, of class
	 * GenGenotypeData.
	 */
	@Test
	public void testIsOnlyContaingSaveProbabilityGenotypes() {
		assertTrue(genotypeData.isOnlyContaingSaveProbabilityGenotypes());
	}

	@Test
	public void testGetSnpVariantByPos() {
		int pos = 14433624;
		GeneticVariant variant = genotypeData.getSnpVariantByPos("22", pos);
		assertNotNull(variant);
		assertEquals(variant.getStartPos(), pos);


		float[][] expResult = new float[9][3];

		expResult[0] = new float[]{0, 0, 1};
		expResult[1] = new float[]{0, 1, 0};
		expResult[2] = new float[]{0, 0, 1};
		expResult[3] = new float[]{0, 0, 1};
		expResult[4] = new float[]{0, 0, 1};
		expResult[5] = new float[]{0, 1, 0};
		expResult[6] = new float[]{0, 0, 1};
		expResult[7] = new float[]{0, 0, 1};
		expResult[8] = new float[]{0, 1, 0};

		assertEquals(variant.getSampleGenotypeProbilities(), expResult, 0.001f, "Probs not identical");

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
		assertEquals(variant.getVariantAlleles(), Alleles.createBasedOnChars('A', 'G'));

		byte[] expectedCalledDosage = new byte[]{0, 1, 0, 0, 0, 1, 0, 0, 1};

		assertEquals(variant.getSampleCalledDosages(), expectedCalledDosage);

	}

	@Test
	public void testGetSnpVariantByPos2() {
		int pos = 14432918;
		GeneticVariant variant = genotypeData.getSnpVariantByPos("22", pos);
		assertNotNull(variant);
		assertEquals(variant.getStartPos(), pos);

		assertEquals(variant.getVariantAlleles(), Alleles.createAlleles(Allele.C));


		float[][] expResult = new float[9][3];

		expResult[0] = new float[]{1, 0, 0};
		expResult[1] = new float[]{1, 0, 0};
		expResult[2] = new float[]{1, 0, 0};
		expResult[3] = new float[]{1, 0, 0};
		expResult[4] = new float[]{1, 0, 0};
		expResult[5] = new float[]{1, 0, 0};
		expResult[6] = new float[]{1, 0, 0};
		expResult[7] = new float[]{1, 0, 0};
		expResult[8] = new float[]{0, 0, 0};

		assertEquals(variant.getSampleGenotypeProbilities(), expResult, 0.001f, "Probs not identical");

		ArrayList<Alleles> expectedSampleAlleles = new ArrayList<Alleles>(9);
		expectedSampleAlleles.add(Alleles.createBasedOnChars('C', 'C'));
		expectedSampleAlleles.add(Alleles.createBasedOnChars('C', 'C'));
		expectedSampleAlleles.add(Alleles.createBasedOnChars('C', 'C'));
		expectedSampleAlleles.add(Alleles.createBasedOnChars('C', 'C'));
		expectedSampleAlleles.add(Alleles.createBasedOnChars('C', 'C'));
		expectedSampleAlleles.add(Alleles.createBasedOnChars('C', 'C'));
		expectedSampleAlleles.add(Alleles.createBasedOnChars('C', 'C'));
		expectedSampleAlleles.add(Alleles.createBasedOnChars('C', 'C'));
		expectedSampleAlleles.add(Alleles.createBasedOnChars('0', '0'));

		assertEquals(variant.getSampleVariants(), expectedSampleAlleles);
		
		float[] expectedDosage = new float[]{2f, 2f, 2f, 2f, 2f, 2f, 2f, 2f, -1f};

		assertEquals(variant.getSampleDosages(), expectedDosage, 0.0001f, "Dosage not identical");
		

		byte[] expectedCalledDosage = new byte[]{2, 2, 2, 2, 2, 2, 2, 2, -1};


		assertEquals(variant.getSampleCalledDosages(), expectedCalledDosage);

	}

	@Test
	public void testGetSnpVariantByPos3() {
		int pos = 1;
		GeneticVariant variant = genotypeData.getSnpVariantByPos("1", pos);
		assertNotNull(variant);
		assertEquals(variant.getStartPos(), pos);

		assertEquals(variant.getVariantAlleles(), Alleles.createAlleles(Allele.A, Allele.T));


		float[][] expResult = new float[9][3];

		expResult[0] = new float[]{0.5f, 0.25f, 0.25f};
		expResult[1] = new float[]{0.25f, 0.5f, 0.25f};
		expResult[2] = new float[]{0.25f, 0.25f, 0.5f};
		expResult[3] = new float[]{0.1f, 0.1f, 0.1f};
		expResult[4] = new float[]{0.3f, 0.3f, 0.3f};
		expResult[5] = new float[]{0.4f, 0.3f, 0.3f};
		expResult[6] = new float[]{0.2f, 0.3f, 0.5f};
		expResult[7] = new float[]{0, 0.1f, 1f};
		expResult[8] = new float[]{0.01f, 0.001f, 0.989f};

		assertEquals(variant.getSampleGenotypeProbilities(), expResult, 0.001f, "Probs not identical");

		ArrayList<Alleles> expectedSampleAlleles = new ArrayList<Alleles>(9);
		expectedSampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		expectedSampleAlleles.add(Alleles.createBasedOnChars('A', 'T'));
		expectedSampleAlleles.add(Alleles.createBasedOnChars('T', 'T'));
		expectedSampleAlleles.add(Alleles.createBasedOnChars('0', '0'));
		expectedSampleAlleles.add(Alleles.createBasedOnChars('0', '0'));
		expectedSampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		expectedSampleAlleles.add(Alleles.createBasedOnChars('T', 'T'));
		expectedSampleAlleles.add(Alleles.createBasedOnChars('T', 'T'));
		expectedSampleAlleles.add(Alleles.createBasedOnChars('T', 'T'));

		assertEquals(variant.getSampleVariants(), expectedSampleAlleles);

		float[] expectedDosage = new float[]{1.25f, 1f, 0.75f, -1f, -1f, 1.1f, 0.7f, 0.1f, 0.021f};

		assertEquals(variant.getSampleDosages(), expectedDosage, 0.0001f, "Dosage not identical");

		byte[] expectedCalledDosage = new byte[]{2, 1, 0, -1, -1, 2, 0, 0, 0};

		assertEquals(variant.getSampleCalledDosages(), expectedCalledDosage);

	}
	
}