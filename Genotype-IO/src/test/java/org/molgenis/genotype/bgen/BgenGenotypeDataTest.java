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
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;


import org.apache.commons.io.FilenameUtils;
import org.junit.rules.TemporaryFolder;
import org.molgenis.genotype.ResourceTest;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.oxford.GenGenotypeData;
import org.molgenis.genotype.oxford.HapsGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.ReadOnlyGeneticVariantBgen;
import org.molgenis.genotype.vcf.VcfGenotypeData;
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
	private File exampleSampleFile = getTestResourceFile("/bgenExamples/genFiles/example.sample");
	private File exampleGenFile = getTestResourceFile("/bgenExamples/genFiles/example.gen");
	private File complexVcfFile = getTestResourceFile("/bgenExamples/complex.vcf.gz");
	private List<File> complexFiles = Arrays.asList(
			getTestResourceFile("/bgenExamples/complex.3bits.bgen"),
			getTestResourceFile("/bgenExamples/complex.5bits.bgen"),
			getTestResourceFile("/bgenExamples/complex.15bits.bgen"),
			getTestResourceFile("/bgenExamples/complex.24bits.bgen"),
			getTestResourceFile("/bgenExamples/complex.31bits.bgen"),
			getTestResourceFile("/bgenExamples/complex.24bits.bgen"),
			getTestResourceFile("/bgenExamples/complex.bgen")
	);
 	private List<File> exampleFiles = Arrays.asList(
			getTestResourceFile("/bgenExamples/example.1bits.bgen"),
			getTestResourceFile("/bgenExamples/example.2bits.bgen"),
			getTestResourceFile("/bgenExamples/example.8bits.bgen"),
			getTestResourceFile("/bgenExamples/example.16bits.bgen"),
			getTestResourceFile("/bgenExamples/example.16bits.zstd.bgen"),
			getTestResourceFile("/bgenExamples/example.25bits.bgen"),
			getTestResourceFile("/bgenExamples/example.32bits.bgen")
//			getTestResourceFile("/bgenExamples/example.v11.bgen")
	);

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
	public void bgenGenotypeDataTest() throws IOException {

		// Load the example .gen file corresponding to the BGEN example files.
		GenGenotypeData genGenotypeData = new GenGenotypeData(exampleGenFile, exampleSampleFile);

		// Loop through the bgen example files.
		for (File origBgenFile : exampleFiles) {

			// Try to extract the bit representation of the example file my matching a regex with a group for the
			// bit representation number.
			Matcher matcher = Pattern.compile("example\\.(\\d+)bits\\.(zstd\\.)?bgen")
					.matcher(origBgenFile.getName());
			// Get the bit representation.
			int bitRepresentation = matcher.matches() ? Integer.parseInt(matcher.group(1)) : 0;

			// Set the maximum error from the bitrepresentation as shown within the BGEN specification.
			double maximumError = 1 / (Math.pow(2, bitRepresentation) - 1);
			// Replace the maximum error by 1*10^4 if the original maximum error was lower, as the precision of
			// the probabilities in the .gen file is often not greater than 4 decimals.
			maximumError = Math.max(maximumError, 1e-4);

			// Load the bgen file from a temporary folder
			Path bgenFile = Paths.get(folder.getRoot().toString(), origBgenFile.getName());
			Files.copy(origBgenFile.toPath(), bgenFile);
			bgenGenotypeData = new BgenGenotypeData(bgenFile.toFile());

			// Test the equality of sample names and sequence names
			assertEquals(bgenGenotypeData.getSampleNames(), genGenotypeData.getSampleNames());
			assertEquals(bgenGenotypeData.getSeqNames(), genGenotypeData.getSeqNames());

			Iterator<GeneticVariant> bgenIterator = bgenGenotypeData.iterator();
			// Loop through variants and check their similarity.
			for (GeneticVariant genVariant : genGenotypeData) {

				// Check if the bgen iterator also has a next variant.
				assertTrue(bgenIterator.hasNext(), "bgenIterator is emptied while Oxford genIterator is not.");
				GeneticVariant bgenVariant = bgenIterator.next();

				// Check if these are equal:
				assertEquals(bgenVariant.getPrimaryVariantId(), genVariant.getPrimaryVariantId());
				assertEquals(bgenVariant.getVariantAlleles(), genVariant.getVariantAlleles());
				assertEquals(bgenVariant.getAlternativeAlleles(), genVariant.getAlternativeAlleles());
				assertEquals(bgenVariant.getRefAllele(), genVariant.getRefAllele());
				assertEquals(bgenVariant.getSequenceName(), genVariant.getSequenceName());
				assertEquals(bgenVariant.getStartPos(), genVariant.getStartPos());
				assertEquals(bgenVariant.getAlleleCount(), genVariant.getAlleleCount());

				// Check if the probabilities are similar enough to the probabilities within the .gen file
				float[][] bgenProbabilities = bgenVariant.getSampleGenotypeProbilities();

				// Since the variant RSID_10 is effectively equal to RSID_100 within the gen file (same chr and bp),
				// the genotype probabilities returned by the gen file for RSID_10 correspond to RSID_100 instead.
				// Therefore we check some probabilities for this variant manually.
				if (genVariant.getPrimaryVariantId().equals("RSID_10")) {
					float[][] genProbabilities = {
							{0.0145569f, 0.960388f, 0.0250549f},
							{0.0101318f, 0.989807f, 6.10352e-05f,}};
					bgenProbabilities = Arrays.copyOfRange(bgenProbabilities, 0, 2);
					assertProbabilityEquality(bgenProbabilities, genProbabilities, maximumError);
				} else {
					// If the variant is not equal to RSID_10, just use the probabilities returned by the
					// genVariant.
					float[][] genProbabilities = genVariant.getSampleGenotypeProbilities();
					assertProbabilityEquality(bgenProbabilities, genProbabilities, maximumError);
				}
			}
		}
	}

	private static void assertProbabilityEquality(float[][] actual, float[][] expected, double maximumError) {
		assertEquals(actual.length, expected.length);
		for(int i = 0; i < actual.length; i++){
			assertEquals(actual[i].length, expected[i].length,
					String.format("number of probabilities not equal for sample %d", i));
			for(int y = 0; y < actual[i].length; y++){
				assertEquals(actual[i][y], expected[i][y],
						maximumError, String.format("prob %d in sample %d", y, i));
			}
		}
	}

	private static void assertProbabilityEquality(double[][] actual, double[][] expected) {
		assertEquals(actual.length, expected.length);
		for(int i = 0; i < actual.length; i++){
			assertEquals(actual[i].length, expected[i].length,
					String.format("number of probabilities not equal for sample %d", i));
			for(int y = 0; y < actual[i].length; y++){
				assertEquals(actual[i][y], expected[i][y], String.format("prob %d in sample %d", y, i));
			}
		}
	}

	@Test
	public void complexBgenGenotypeDataTest() throws IOException {

		// Load the example VCF file corresponding to the complex[.<number of bits>bits].bgen files
		VcfGenotypeData vcfGenotypeData = new VcfGenotypeData(
				complexVcfFile, 0, 0.4f);

		// Loop through the complex bgen files and load the corresponding copy from a temporary folder
		for (File origBgenFile : complexFiles) {
			// Copy the BGEN file and load this copy
			Path bgenFile = Paths.get(folder.getRoot().toString(), origBgenFile.getName());
			Files.copy(origBgenFile.toPath(), bgenFile);
			bgenGenotypeData = new BgenGenotypeData(bgenFile.toFile());

			// Test the equality of sample names and sequence names
			assertEquals(bgenGenotypeData.getSampleNames(), vcfGenotypeData.getSampleNames());
			assertEquals(bgenGenotypeData.getSeqNames(), vcfGenotypeData.getSeqNames());

			Iterator<GeneticVariant> bgenIterator = bgenGenotypeData.iterator();
			// Loop through variants and check their similarity.
			for (GeneticVariant vcfVariant : vcfGenotypeData) {
				// Check if the bgen iterator also has a next variant.
				assertTrue(bgenIterator.hasNext(), "bgenIterator is emptied while Oxford genIterator is not.");
				GeneticVariant bgenVariant = bgenIterator.next();

				// Check if these variants are equal:
				assertEquals(bgenVariant.getPrimaryVariantId(), vcfVariant.getPrimaryVariantId());
				assertEquals(bgenVariant.getSequenceName(), vcfVariant.getSequenceName());
				assertEquals(bgenVariant.getStartPos(), vcfVariant.getStartPos());
				assertEquals(bgenVariant.getAlleleCount(), vcfVariant.getAlleleCount());

				// Check the equality of probabilities.
				double[][] bgenProbabilities = bgenVariant.getSampleGenotypeProbabilitiesBgen();
				System.out.println(Arrays.deepToString(bgenProbabilities));
				double[][] vcfProbabilities = vcfVariant.getSampleGenotypeProbabilitiesBgen();
				System.out.println(Arrays.deepToString(vcfProbabilities));
				assertProbabilityEquality(bgenProbabilities, vcfProbabilities);
			}

		}
	}

//	@Test
//	public void testReader_1_2() throws URISyntaxException, IOException {
//		//File bgenFile = getTestResourceFile("/bgenExamples/complex.bgen");
//		//genotypeData = new BgenGenotypeData(bgenFile, null);
//		File bgenFile = getTestResourceFile("/bgenExamples/example.25bits.bgen");
//		Path target = Paths.get(folder.getRoot().toString(), bgenFile.getName());
//		Files.copy(bgenFile.toPath(), target);
//		bgenGenotypeData = new BgenGenotypeData(target.toFile());
//
//		int counter = 0;
//		for (GeneticVariant variant : bgenGenotypeData) {
//			System.out.println("variant.getPrimaryVariantId() = " + variant.getPrimaryVariantId());
//			float[][] sampleGenotypeProbilities = variant.getSampleGenotypeProbilities();
//			System.out.println("sampleGenotypeProbilities = " + Arrays.deepToString(sampleGenotypeProbilities));
//		}
//	}

//	@Test
//	public void testReader_1_3() throws URISyntaxException, IOException {
//		File bgenFile = getTestResourceFile("/bgenExamples/example.16bits.zstd.bgen");
//		Path target = Paths.get(folder.getRoot().toString(), bgenFile.getName());
//		Files.copy(bgenFile.toPath(), target);
//		bgenGenotypeData = new BgenGenotypeData(target.toFile());
//	}

	@Test
	public void testReaderHaplotypes() throws URISyntaxException, IOException {
		// Get the bgen input file to test with
		File bgenFile = getTestResourceFile("/bgenExamples/haplotypes.bgen");
		// Get the corresponding haplotype & sample files
		File hapsFile = getTestResourceFile("/bgenExamples/genFiles/haplotypes.haps");
		File oxfordSampleFile = getTestResourceFile("/bgenExamples/genFiles/haplotypes.sample");

		// Copy the input file
		Path target = Paths.get(folder.getRoot().toString(), bgenFile.getName());
		Files.copy(bgenFile.toPath(), target);
		// Load the haps and bgen data
		HapsGenotypeData hapsGenotypeData = new HapsGenotypeData(hapsFile, oxfordSampleFile);
		bgenGenotypeData = new BgenGenotypeData(target.toFile());

		// Test the equality of sample names and sequence names
		assertEquals(bgenGenotypeData.getSampleNames(), hapsGenotypeData.getSampleNames());
		assertEquals(bgenGenotypeData.getSeqNames(), hapsGenotypeData.getSeqNames());

		Iterator<GeneticVariant> bgenIterator = bgenGenotypeData.iterator();
		// Loop trough the variants of both the haps data and the bgen data.
		for (GeneticVariant hapsVariant : hapsGenotypeData) {
			// Check if the number of variants is equal and get the next variant from the bgen iterator.
			assertTrue(bgenIterator.hasNext(), "bgenIterator is emptied while Oxford hapsIterator is not.");
			GeneticVariant bgenVariant = bgenIterator.next();

			// Check for equality between the data objects
			assertEquals(bgenVariant.getPrimaryVariantId(), hapsVariant.getPrimaryVariantId());
			assertEquals(bgenVariant.getSequenceName(), hapsVariant.getSequenceName());
			assertEquals(bgenVariant.getStartPos(), hapsVariant.getStartPos());
			assertEquals(bgenVariant.getAlleleCount(), hapsVariant.getAlleleCount());

			// Compare probabilities
			float[][] bgenProbabilities = bgenVariant.getSampleGenotypeProbilities();
			float[][] hapsProbabilities = hapsVariant.getSampleGenotypeProbilities();
			assertProbabilityEquality(bgenProbabilities, hapsProbabilities, 0);
		}
	}
}
