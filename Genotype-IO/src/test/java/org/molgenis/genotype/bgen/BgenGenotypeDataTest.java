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
import java.util.Arrays;
import java.util.List;


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
			getTestResourceFile("/bgenExamples/example.16bits.zstd.bgen"),
//			getTestResourceFile("/bgenExamples/example.25bits.bgen"),
			getTestResourceFile("/bgenExamples/example.32bits.bgen"));

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
	public void bgenGenotypeDataTest() throws URISyntaxException, IOException {
		for (File origBgenFile : exampleFiles) {
			Path bgenFile = Paths.get(folder.getRoot().toString(), origBgenFile.getName());
			Files.copy(origBgenFile.toPath(), bgenFile);
			bgenGenotypeData = new BgenGenotypeData(bgenFile.toFile());

			Path exampleFiles = origBgenFile.toPath().resolveSibling("genFiles");

			Path genSampleFile = exampleFiles.resolve(
					bgenFile.getFileName().toString().replace(".bgen", ".sample"));
			Path genFile = exampleFiles.resolve(
					bgenFile.getFileName().toString().replace(".bgen", ".gen"));

			GenGenotypeData genGenotypeData = new GenGenotypeData(genFile.toFile(), genSampleFile.toFile());

			assertEquals(bgenGenotypeData.getSamples(), genGenotypeData.getSamples());

//			for (GeneticVariant variant : bgenGenotypeData) {
//				System.out.println("variant.getPrimaryVariantId() = " + variant.getPrimaryVariantId());
//				float[][] sampleGenotypeProbilities = variant.getSampleGenotypeProbilities();
//				System.out.println("sampleGenotypeProbilities = " + Arrays.deepToString(sampleGenotypeProbilities));
//			}

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
			writer.write(String.format("ID_1 ID_2%n"));
			writer.write(String.format("0 0%n"));
			for (Sample sample : bgenGenotypeData.getSamples()) {
				writer.write(String.format("%s %s%n", sample.getId(), sample.getId()));
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
