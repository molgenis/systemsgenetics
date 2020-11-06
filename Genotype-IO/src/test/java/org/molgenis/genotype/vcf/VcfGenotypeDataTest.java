package org.molgenis.genotype.vcf;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.nio.file.Paths;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;

import org.apache.commons.io.FilenameUtils;
import org.molgenis.genotype.*;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.VcfAnnotation;
import org.molgenis.genotype.bgen.BgenGenotypeData;
import org.molgenis.genotype.bgen.BgenGenotypeWriter;
import org.molgenis.genotype.util.Utils;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variantFilter.VariantIdIncludeFilter;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

import com.google.common.collect.Lists;

import static org.testng.Assert.*;

public class VcfGenotypeDataTest extends ResourceTest
{
	private VcfGenotypeData genotypeData;
	private VcfGenotypeData complexGenotypeData;
	private File folder;
	private File complexVcfFile;

	@BeforeClass
	public void beforeClass() throws IOException, URISyntaxException
	{
		genotypeData = new VcfGenotypeData(getTestVcfGz(), getTestVcfGzTbi(), 0.8);
		complexVcfFile = getTestResourceFile("/bgenExamples/complex.vcf.gz");
		File complexVcfIndexFile = getTestResourceFile("/bgenExamples/complex.vcf.gz.tbi");
		complexGenotypeData = new VcfGenotypeData(
				complexVcfFile,
				complexVcfIndexFile, 0.0);
	}

	@BeforeTest
	public void setUpMethod() {
		File tmpDir = new File(System.getProperty("java.io.tmpdir"));

		DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
		Date date = new Date();

		folder = new File(tmpDir, "VcfGenotypeDataTest_" + dateFormat.format(date));

		Runtime.getRuntime().addShutdownHook(new Thread(() -> {
			System.out.println("Removing tmp dir and files");
			File[] files = folder.listFiles();
			if (files != null) {
				for (File file : files) {
					System.out.println(" - Deleting: " + file.getAbsolutePath());
					assert file.delete();
				}
			}
			System.out.println(" - Deleting: " + folder.getAbsolutePath());
			assert folder.delete();
		}));

		assert folder.mkdir();

		System.out.println("Temp folder with output of this test: " + folder.getAbsolutePath());
	}

	@Test
	public void getSeqNames()
	{
		List<String> seqNames = genotypeData.getSeqNames();
		assertNotNull(seqNames);
		assertEquals(seqNames.size(), 3);
		assertEquals(seqNames.get(0), "1");
		assertEquals(seqNames.get(1), "2");
		assertEquals(seqNames.get(2), "3");
	}

	@Test
	public void getVariantAnnotationsMap()
	{
		Map<String, Annotation> annotations = genotypeData.getVariantAnnotationsMap();
		assertNotNull(annotations);
		assertEquals(annotations.size(), 21);

		Annotation annotation = annotations.get("NS");
		assertNotNull(annotation);
		assertEquals(annotation.getId(), "NS");
		assertTrue(annotation instanceof VcfAnnotation);
	}

	@Test
	public void testGetSequences()
	{
		List<Sequence> sequences = genotypeData.getSequences();
		assertNotNull(sequences);
		assertEquals(sequences.size(), 3);
	}

	@Test
	public void testGetSequenceByName() throws IOException
	{
		Sequence sequence = genotypeData.getSequenceByName("2");
		assertNotNull(sequence);
		assertEquals(sequence.getName(), "2");

		List<GeneticVariant> variants = Utils.iteratorToList(genotypeData.getSequenceGeneticVariants(sequence.getName()).iterator());

		assertNotNull(variants);
		assertEquals(variants.size(), 1);
		GeneticVariant variant = variants.get(0);
		assertEquals(variant.getPrimaryVariantId(), "rs4908464");
		assertEquals(variant.getStartPos(), 7569187);
		assertEquals(variant.getRefAllele().getAlleleAsString(), "G");
		assertEquals(variant.getAlternativeAlleles(), Alleles.createAlleles(Allele.C));
		assertEquals(variant.getSequenceName(), "2");
		assertEquals(variant.isSnp(), true);

		List<String> alleles = variant.getVariantAlleles().getAllelesAsString();
		assertNotNull(alleles);
		assertEquals(alleles.size(), 2);
		assertEquals(alleles.get(0), "G");
		assertEquals(alleles.get(1), "C");

		List<Alleles> sampleVariants = variant.getSampleVariants();
		assertNotNull(sampleVariants);
		assertEquals(sampleVariants.size(), 1);
		assertNotNull(sampleVariants.get(0).getAlleles());
		assertEquals(sampleVariants.get(0).getAlleles().size(), 2);
		assertEquals(sampleVariants.get(0).getAlleles().get(0).getAlleleAsString(), "C");
		assertEquals(sampleVariants.get(0).getAlleles().get(1).getAlleleAsString(), "C");
	}

	@Test
	public void testGetSequenceByNameNotExisting()
	{
		assertNull(genotypeData.getSequenceByName("Bogus"));
	}

	@Test
	public void testGetVariant()
	{
		List<GeneticVariant> variants = Lists.newArrayList(genotypeData.getVariantsByPos("1", 565286));
		assertNotNull(variants);
		assertEquals(variants.size(), 1);
		assertEquals(variants.get(0).getPrimaryVariantId(), "rs1578391");

		assertNotNull(genotypeData.getVariantsByPos("bogus", 8));
		assertFalse(genotypeData.getVariantsByPos("bogus", 8).iterator().hasNext());
	}

	@Test
	public void testgetSequenceGeneticVariants() throws IOException
	{
		List<GeneticVariant> variants = Utils.iteratorToList(genotypeData.getSequenceGeneticVariants("1").iterator());
		assertEquals(variants.size(), 6);
	}

	@Test
	public void testSnpVariants()
	{
		GeneticVariant snpGeneticVariant = genotypeData.getSnpVariantByPos("1", 3172273);
		assertNotNull(snpGeneticVariant);
		assertEquals(snpGeneticVariant.isSnp(), true);
	}

	@Test
	public void testGeneticVariantAnnotations()
	{
		// NS=1;DP=4;AF=1.000;ANNOT=INT;GI=PRDM16;TI=NM_022114.3;PI=NP_071397.3
		List<GeneticVariant> variants = Lists.newArrayList(genotypeData.getVariantsByPos("1", 3171929));
		assertNotNull(variants);
		assertEquals(variants.size(), 1);

		GeneticVariant variant = variants.get(0);
		assertNotNull(variant.getAnnotationValues());
		assertEquals(variant.getAnnotationValues().size(), 9);//9 because of filter and qual field

		Object annotationValue = variant.getAnnotationValues().get("NS");
		assertNotNull(annotationValue);
		assertEquals(annotationValue, Integer.valueOf(1));

		annotationValue = variant.getAnnotationValues().get("AF");
		assertNotNull(annotationValue);
		assertTrue(annotationValue instanceof List);
		@SuppressWarnings("unchecked")
		List<Float> floats = (List<Float>) annotationValue;
		assertEquals(floats.size(), 1);
		assertEquals(floats.get(0).floatValue(), 1.0, 0.001);

		annotationValue = variant.getAnnotationValues().get("ANNOT");
		assertNotNull(annotationValue);
		assertTrue(annotationValue instanceof List);
		@SuppressWarnings("unchecked")
		List<String> strings = (List<String>) annotationValue;
		assertEquals(strings.size(), 1);
		assertEquals(strings.get(0), "INT");
		
		assertEquals(variant.getAnnotationValues().get("VCF_Filter"), "flt");
		assertNull(variant.getAnnotationValues().get("VCF_Qual"));

	}

	@Test
	public void testStopPos() throws IOException, URISyntaxException
	{
		List<GeneticVariant> variants = Lists.newArrayList(genotypeData.getVariantsByPos("1", 565286));
		assertNotNull(variants);
		assertEquals(variants.size(), 1);

		variants = Lists.newArrayList(genotypeData.getVariantsByPos("3", 7569));
		assertNotNull(variants);
		assertEquals(variants.size(), 1);
	}

	@Test
	public void testGetSamplePhasing()
	{
		GeneticVariant variant = genotypeData.getVariantsByPos("1", 565286).iterator().next();
		assertEquals(genotypeData.getSamplePhasing(variant), Arrays.asList(false));
	}
	
	@Test
	public void getVariantIdMap(){
		HashMap<String, GeneticVariant> variantMap = genotypeData.getVariantIdMap();
		
		assertEquals(variantMap.size(), 7);
		
		assertTrue(variantMap.containsKey("rs35434908"));
		
		assertEquals(variantMap.get("rs35434908").getVariantAlleles().get(1), Allele.create("GTTTCA"));
		
		variantMap = genotypeData.getVariantIdMap(new VariantIdIncludeFilter("rs35434908", "rs1578391"));

		assertEquals(variantMap.size(), 2);
		
		assertTrue(variantMap.containsKey("rs35434908"));
		assertTrue(variantMap.containsKey("rs1578391"));
		
	}


	@Test
	public void complexVcfGenotypeDataTest() throws IOException {

		// Initialize the expected stuff...
		String[] expectedSampleNames = {"sample_0", "sample_1", "sample_2", "sample_3"};
		List<String> expectedVariantIds = Arrays.asList("V1", "V2", "V3", "M4", "M5", "M6", "M7", "M8", "M9", "M10");
		List<Integer> expectedVariantPositions = Arrays.asList(1, 2, 3, 4, 5, 7, 7, 8, 9, 10);
		List<Alleles> alleles = Arrays.asList(
				Alleles.createAlleles(Allele.A, Allele.G),
				Alleles.createAlleles(Allele.A, Allele.G),
				Alleles.createAlleles(Allele.A, Allele.G),
				Alleles.createAlleles(Allele.A, Allele.G, Allele.T),
				Alleles.createAlleles(Allele.A, Allele.G),
				Alleles.createAlleles(Allele.A, Allele.G, Allele.create("GT"), Allele.create("GTT")),
				Alleles.createAlleles(Allele.A, Allele.G, Allele.create("GT"), Allele.create("GTT"), Allele.create("GTTT"), Allele.create("GTTTT")),
				Alleles.createAlleles(Allele.A, Allele.G, Allele.create("GT"), Allele.create("GTT"), Allele.create("GTTT"), Allele.create("GTTTT"), Allele.create("GTTTTT")),
				Alleles.createAlleles(Allele.A, Allele.G, Allele.create("GT"), Allele.create("GTT"), Allele.create("GTTT"), Allele.create("GTTTT"), Allele.create("GTTTTT"), Allele.create("GTTTTTT")),
				Alleles.createAlleles(Allele.A, Allele.G));

		// Initialize interesting probability arrays to check against
		List<double[][]> expectedComplexProbabilities = Arrays.asList(
				new double[][]{{1, 0}, {1, 0, 0}, {1, 0, 0}, {0, 1, 0}},
				new double[][]{{1, 0}, {1, 0}, {1, 0}, {0, 1}},
				new double[][]{{1, 0}, {0, 0, 1}, {1, 0, 0}, {0, 1, 0}},
				new double[][]{{1, 0, 0}, {1, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0}},
				new double[][]{{1, 0}, {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 1, 0}},
				null,
				null,
				null,
				new double[][]{{1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
				new double[][]{{1, 0, 0, 0, 0}, {0, 1, 0, 0, 0}, {0, 0, 1, 0, 0}, {0, 0, 0, 1, 0}});

		List<float[][]> expectedProbabilities = Arrays.asList(
				new float[][]{{0, 0, 0}, {1, 0, 0}, {1, 0, 0}, {0, 1, 0}},
				null,
				null,
				new float[][]{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
				new float[][]{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 1, 0}});

		List<double[][][]> expectedPhasedProbabilities = Arrays.asList(
				null,
				new double[][][]{{{1, 0}}, {{1, 0}}, {{1, 0}}, {{0, 1}}},
				new double[][][]{{{1, 0}}, {{0, 1}, {0, 1}}, {{1, 0}, {1, 0}}, {{1, 0}, {0, 1}}},
				null,
				new double[][][]{{{1, 0}}, {{1, 0}, {1, 0}, {1, 0}}, {{1, 0}, {1, 0}, {0, 1}}, {{1, 0}, {0, 1}}},
				new double[][][]{{{1, 0, 0, 0}}, {{0, 1, 0, 0}}, {{0, 0, 1, 0}}, {{0, 0, 0, 1}}});

		testComplexGenotypeDataCorrespondence(expectedSampleNames, expectedVariantIds, expectedVariantPositions, alleles,
				expectedComplexProbabilities, expectedProbabilities, expectedPhasedProbabilities, complexGenotypeData);
		int variantIndex;

		variantIndex = 0;
		for (GeneticVariant variant : complexGenotypeData) {
			complexGenotypeData.setPreferredGenotypeField("GT");

			// Check the equality of probabilities.
			// First check if the bgenProbabilities are according to the expected stuff
			if (Arrays.asList(0, 1, 2, 3, 4, 8, 9).contains(variantIndex)) {
				double[][] complexProbabilities = variant.getSampleGenotypeProbabilitiesComplex();
				assertEquals(complexProbabilities, expectedComplexProbabilities.get(variantIndex));
			}
			variantIndex++;
		}

		complexGenotypeData.setPreferredGenotypeField(null);

		// Write and read again
		BgenGenotypeWriter bgenGenotypeWriter = new BgenGenotypeWriter(complexGenotypeData);

		File tempFile = Paths.get(folder.toString(),
				FilenameUtils.removeExtension(complexVcfFile.getName()) + ".temp.bgen").toFile();
		File tempSampleFile = Paths.get(folder.toString(),
				FilenameUtils.removeExtension(complexVcfFile.getName()) + ".temp.sample").toFile();
		bgenGenotypeWriter.write(tempFile, tempSampleFile);
		BgenGenotypeData bgenGenotypeData = new BgenGenotypeData(tempFile, tempSampleFile);

		testComplexGenotypeDataCorrespondence(expectedSampleNames,
				expectedVariantIds, expectedVariantPositions,
				alleles, expectedComplexProbabilities,
				expectedProbabilities, expectedPhasedProbabilities, bgenGenotypeData);
	}

	private void testComplexGenotypeDataCorrespondence(String[] expectedSampleNames, List<String> expectedVariantIds,
													   List<Integer> expectedVariantPositions, List<Alleles> alleles,
													   List<double[][]> expectedComplexProbabilities,
													   List<float[][]> expectedProbabilities,
													   List<double[][][]> expectedPhasedProbabilities,
													   RandomAccessGenotypeData complexGenotypeData) {

		// Test the equality of sample names and sequence names
		assertEquals(this.complexGenotypeData.getSampleNames(), expectedSampleNames);
		assertEquals(this.complexGenotypeData.getSeqNames(), Collections.singletonList("01"));

		assertEquals(this.complexGenotypeData.getVariantIdMap().keySet(), new HashSet<>(expectedVariantIds));
		// Loop through variants and check their similarity.
		int variantIndex = 0;
		for (GeneticVariant variant : complexGenotypeData) {
			// Check if these variants are equal:
			assertEquals(variant.getPrimaryVariantId(), expectedVariantIds.get(variantIndex));
			assertEquals(variant.getStartPos(), expectedVariantPositions.get(variantIndex).intValue());
			assertEquals(variant.getVariantAlleles(), alleles.get(variantIndex));

			// V2 has an alternative id, the other variants don't. Check this.
			if (variant.getPrimaryVariantId().equals("V2")) {
				assertEquals(variant.getAlternativeVariantIds().get(0), "V2.1");
			} else {
				assertEquals(variant.getAlternativeVariantIds().size(), 0);
			}

			// Check the equality of probabilities.
			// First check if the bgenProbabilities are according to the expected stuff
			if (Arrays.asList(0, 1, 2, 3, 4, 8, 9).contains(variantIndex)) {
				double[][] complexProbabilities = variant.getSampleGenotypeProbabilitiesComplex();
				assertEquals(complexProbabilities, expectedComplexProbabilities.get(variantIndex));
			}

			// Check if the regular probabilities are according to the expected stuff
			// (a lot should be coded as being missing)
			if (Arrays.asList(0, 3, 4).contains(variantIndex)) {
				float[][] probabilities = variant.getSampleGenotypeProbilities();
				assertEquals(
						probabilities,
						expectedProbabilities.get(variantIndex));
			}

			// Check if the phased probabilities are correct as well
			if (Arrays.asList(2, 4).contains(variantIndex)) {
				// First we have to check if the variant is phased
				assertTrue(variant.hasPhasedProbabilities());
				// Only then get the phased probabilities
				double[][][] phasedBgenProbabilities = variant.getSampleGenotypeProbabilitiesPhased();
				assertTrue(Arrays.deepEquals(phasedBgenProbabilities, expectedPhasedProbabilities.get(variantIndex)));
			} else {
				// These are not phased, which we have to test for
				assertFalse(variant.getSamplePhasing().contains(true));
				assertFalse(variant.hasPhasedGenotypes());
				assertFalse(variant.hasPhasedProbabilities());
				// Calling the method then should generate an exception
				try {
					variant.getSampleGenotypeProbabilitiesPhased();
					fail("variant.getSampleGenotypeProbabilitiesPhased() did not raise a "
							+ "GenotypeDataException while phased data was not available");
				} catch (GenotypeDataException e) {
					assertEquals(e.getMessage(), "Phased data not available");
				}
			}
			variantIndex++;
		}
	}

}
