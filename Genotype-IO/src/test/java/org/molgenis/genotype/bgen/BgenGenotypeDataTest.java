/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.bgen;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;


import org.junit.rules.TemporaryFolder;
import org.molgenis.genotype.ResourceTest;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.oxford.GenGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.ReadOnlyGeneticVariantBgen;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import static org.testng.Assert.*;

/**
 *
 * @author Patrick Deelen
 */
public class BgenGenotypeDataTest extends ResourceTest {

	private BgenGenotypeData bgenGenotypeData;
	private TemporaryFolder folder;
	private List<File> exampleFiles = Arrays.asList(
//			getTestResourceFile("/bgenExamples/example.1bits.bgen"),
//			getTestResourceFile("/bgenExamples/example.2bits.bgen"),
//			getTestResourceFile("/bgenExamples/example.8bits.bgen"),
//			getTestResourceFile("/bgenExamples/example.16bits.bgen"),
			getTestResourceFile("/bgenExamples/example.16bits.zstd.bgen"));
//			getTestResourceFile("/bgenExamples/example.25bits.bgen"),
//			getTestResourceFile("/bgenExamples/example.32bits.bgen"));

	public BgenGenotypeDataTest() throws URISyntaxException {
	}

	@BeforeClass
	public void BeforeClass() throws IOException {
		folder = new TemporaryFolder();
		folder.create();
	}

	@Test
	public void testSomeMethod() {
		// TODO review the generated test code and remove the default call to fail.
	}

	@Test
	public void testBgenGenotypeData() throws URISyntaxException, IOException {
		File bgenFile = getTestResourceFile("/bgenExamples/example.25bits.bgen");
		Path target = Paths.get(folder.getRoot().toString(), bgenFile.getName());
		Files.copy(bgenFile.toPath(), target);
		bgenGenotypeData = new BgenGenotypeData(target.toFile());
		for (GeneticVariant variant : bgenGenotypeData) {
			System.out.printf("%s %s %d %s%n", variant.getPrimaryVariantId(), variant.getSequenceName(), variant.getStartPos(), variant.getVariantAlleles());
		}
	}

	@Test
	public void bgenGenotypeDataTest() throws URISyntaxException, IOException {
		for (File origBgenFile : exampleFiles) {
			Path bgenFile = Paths.get(folder.getRoot().toString(), origBgenFile.getName());
			Files.copy(origBgenFile.toPath(), bgenFile);
			bgenGenotypeData = new BgenGenotypeData(bgenFile.toFile());

			Path exampleFiles = origBgenFile.toPath().resolveSibling("genFiles");

			Path genSampleFile = exampleFiles.resolve(
					bgenFile.getFileName().toString().replace(".bgen", ".sample"));
			Path genFile = exampleFiles.resolve(
					bgenFile.getFileName().toString().replace(".bgen", ".a.gen"));

			GenGenotypeData genGenotypeData = new GenGenotypeData(genFile.toFile(), genSampleFile.toFile());

			assertEquals(
					bgenGenotypeData.getSamples().stream().map(Sample::getId).collect(Collectors.toList()),
					genGenotypeData.getSamples().stream().map(Sample::getId).collect(Collectors.toList()));

			Iterator<GeneticVariant> bgenIterator = bgenGenotypeData.iterator();

			ArrayList<GeneticVariant> genVariants = new ArrayList<>();
			for (GeneticVariant genVariant : genGenotypeData) {
				if (!genVariants.contains(genVariant)) {
					genVariants.add(genVariant);
				} else {
					genVariants.remove(genVariant);
				}
			}
			// Loop through variants and check their similarity.
			for (GeneticVariant genVariant : genVariants) {
				// Check if the bgen iterator also has a next variant.
				assertTrue(bgenIterator.hasNext(), "bgenIterator is emptied while Oxford genIterator is not.");
				GeneticVariant bgenVariant = bgenIterator.next();
				// Check if these are equal:
				assertEquals(bgenVariant.getPrimaryVariantId(), genVariant.getPrimaryVariantId());
//				assertEquals(bgenVariant.getPrimaryVariantId(), genVariant.getAllIds());
				assertEquals(bgenVariant.getVariantAlleles(), genVariant.getVariantAlleles());
				assertEquals(bgenVariant.getAlternativeAlleles(), genVariant.getAlternativeAlleles());
				assertEquals(bgenVariant.getRefAllele(), genVariant.getRefAllele());
				assertEquals(bgenVariant.getSequenceName(), genVariant.getSequenceName());
				assertEquals(bgenVariant.getStartPos(), genVariant.getStartPos());
				assertEquals(bgenVariant.getAlleleCount(), genVariant.getAlleleCount());
				float[][] bgenProbabilities = bgenVariant.getSampleGenotypeProbilities();
				float[][] genProbabilities = genVariant.getSampleGenotypeProbilities();
				assertProbabilityEquality(bgenProbabilities, genProbabilities);
			}

		}
	}

	private static void assertProbabilityEquality(float[][] actual, float[][] expected) {
		assertEquals(actual.length, expected.length);
		for(int i = 0; i < actual.length; i++){
			assertEquals(actual[i].length, expected[i].length);
			for(int y = 0; y < actual[i].length; y++){
				assertEquals(actual[i][y], expected[i][y],
						1e-4, String.format("prob %d in sample %d", y, i));
			}
		}
	}

	@Test
	public void writeSampleFiles() throws IOException {
		for (File origBgenFile : exampleFiles) {
			bgenGenotypeData = new BgenGenotypeData(origBgenFile);

			Path bgenFile = origBgenFile.toPath();

			Path exampleFiles = bgenFile.resolveSibling("genFiles");

			Path genSampleFile = exampleFiles.resolve(
					bgenFile.getFileName().toString().replace(".bgen", ".sample"));

			BufferedWriter writer = new BufferedWriter(new FileWriter(genSampleFile.toFile()));
			writer.write(String.format("ID_1 ID_2 missing%n"));
			writer.write(String.format("0 0 0%n"));
			for (Sample sample : bgenGenotypeData.getSamples()) {
				writer.write(String.format("%s %s 0%n", sample.getId(), sample.getId()));
			}
			writer.close();
		}
	}

	@Test
	public void testReader_1_2() throws URISyntaxException, IOException {
		//File bgenFile = getTestResourceFile("/bgenExamples/complex.bgen");
		//genotypeData = new BgenGenotypeData(bgenFile, null);
		File bgenFile = getTestResourceFile("/bgenExamples/example.25bits.bgen");
		Path target = Paths.get(folder.getRoot().toString(), bgenFile.getName());
		Files.copy(bgenFile.toPath(), target);
		bgenGenotypeData = new BgenGenotypeData(target.toFile());

		int counter = 0;
		for (GeneticVariant variant : bgenGenotypeData) {
			System.out.println("variant.getPrimaryVariantId() = " + variant.getPrimaryVariantId());
			float[][] sampleGenotypeProbilities = variant.getSampleGenotypeProbilities();
			System.out.println("sampleGenotypeProbilities = " + Arrays.deepToString(sampleGenotypeProbilities));
		}
	}

	@Test
	public void testReader_1_3() throws URISyntaxException, IOException {
		File bgenFile = getTestResourceFile("/bgenExamples/example.16bits.zstd.bgen");
		Path target = Paths.get(folder.getRoot().toString(), bgenFile.getName());
		Files.copy(bgenFile.toPath(), target);
		bgenGenotypeData = new BgenGenotypeData(target.toFile());
	}

	@Test
	public void testReaderHaplotypes() throws URISyntaxException, IOException {
		File bgenFile = getTestResourceFile("/bgenExamples/haplotypes.bgen");
		Path target = Paths.get(folder.getRoot().toString(), bgenFile.getName());
		Files.copy(bgenFile.toPath(), target);
		bgenGenotypeData = new BgenGenotypeData(target.toFile());

		double[][] trueProbabilitiesBgen = {
				{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
		float[][] trueProbabilities = {
				{1.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f}, {0.0f, 1.0f, 0.0f}, {0.0f, 0.0f, 1.0f}};

		GeneticVariant variant = bgenGenotypeData.iterator().next();
		assertTrue(variant instanceof ReadOnlyGeneticVariantBgen);
		ReadOnlyGeneticVariantBgen bgenVariant = (ReadOnlyGeneticVariantBgen) variant;
		System.out.println("variant = " + bgenVariant.getPrimaryVariantId());
		assertEquals(variant.getPrimaryVariantId(), "RS1");
		double[][] sampleGenotypeProbabilitiesBgen = bgenGenotypeData
				.getSampleGenotypeProbabilitiesBgen(bgenVariant);
		assertTrue(Arrays.deepEquals(trueProbabilitiesBgen, sampleGenotypeProbabilitiesBgen));
		float[][] sampleGenotypeProbabilities = bgenVariant
				.getSampleGenotypeProbilities();
		assertTrue(Arrays.deepEquals(trueProbabilities, sampleGenotypeProbabilities));
	}
}
