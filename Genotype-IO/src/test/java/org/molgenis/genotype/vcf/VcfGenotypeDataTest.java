package org.molgenis.genotype.vcf;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertFalse;
import static org.testng.Assert.assertNotNull;
import static org.testng.Assert.assertNull;
import static org.testng.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;

import org.apache.commons.io.FilenameUtils;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.ResourceTest;
import org.molgenis.genotype.Sequence;
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

		folder = new File(tmpDir, "BgenGenotypeDataTest_" + dateFormat.format(date));

		Runtime.getRuntime().addShutdownHook(new Thread(() -> {
			System.out.println("Removing tmp dir and files");
			for (File file : folder.listFiles()) {
				System.out.println(" - Deleting: " + file.getAbsolutePath());
				file.delete();
			}
			System.out.println(" - Deleting: " + folder.getAbsolutePath());
			folder.delete();
		}));

		folder.mkdir();

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

//		// Write and read again
//		BgenGenotypeWriter bgenGenotypeWriter = new BgenGenotypeWriter(complexGenotypeData);
//
//		File tempFile = Paths.get(folder.toString(),
//				FilenameUtils.removeExtension(complexVcfFile.getName()) + ".temp.bgen").toFile();
//		File tempSampleFile = Paths.get(folder.toString(),
//				FilenameUtils.removeExtension(complexVcfFile.getName()) + ".temp.sample").toFile();
//		bgenGenotypeWriter.write(tempFile, tempSampleFile);
//		BgenGenotypeData bgenGenotypeData = new BgenGenotypeData(tempFile, tempSampleFile);
//		TestComplexBgenGenotypeDataEquality(expectedSampleNames,
//				expectedVariantIds, expectedVariantPositions,
//				alleles, expectedBgenProbabilities,
//				expectedProbabilities, expectedPhasedBgenProbabilities);
	}
	
}
