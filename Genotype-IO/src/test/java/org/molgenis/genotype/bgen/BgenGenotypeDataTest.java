/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.bgen;

import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.ResourceTest;
import org.molgenis.genotype.oxford.GenGenotypeData;
import org.molgenis.genotype.oxford.HapsGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static org.testng.Assert.*;

/**
 * @author Patrick Deelen
 */
public class BgenGenotypeDataTest extends ResourceTest {

    private BgenGenotypeData bgenGenotypeData;
    private File folder;
    private File exampleSampleFile = getTestResourceFile("/bgenExamples/genFiles/example.sample");
    private File exampleComplexSampleFile = getTestResourceFile("/bgenExamples/genFiles/haplotypes.sample");
    private File exampleGenFile = getTestResourceFile("/bgenExamples/genFiles/example.gen");
    private List<File> complexFiles = Arrays.asList(
            getTestResourceFile("/bgenExamples/complex.3bits.bgen"),
            getTestResourceFile("/bgenExamples/complex.5bits.bgen"),
            getTestResourceFile("/bgenExamples/complex.15bits.bgen"),
            getTestResourceFile("/bgenExamples/complex.24bits.bgen"),
            getTestResourceFile("/bgenExamples/complex.31bits.bgen"),
            getTestResourceFile("/bgenExamples/complex.bgen")
    );
    private List<File> exampleFiles = Arrays.asList(
            getTestResourceFile("/bgenExamples/example.1bits.bgen"),
            getTestResourceFile("/bgenExamples/example.2bits.bgen"),
            getTestResourceFile("/bgenExamples/example.8bits.bgen"),
            getTestResourceFile("/bgenExamples/example.16bits.bgen"),
            getTestResourceFile("/bgenExamples/example.16bits.zstd.bgen"),
            getTestResourceFile("/bgenExamples/example.25bits.bgen"),
            getTestResourceFile("/bgenExamples/example.32bits.bgen"),
			getTestResourceFile("/bgenExamples/example.v11.bgen")
    );

    public BgenGenotypeDataTest() throws URISyntaxException {
    }

    @BeforeTest
    public void setUpMethod() {
        File tmpDir = new File(System.getProperty("java.io.tmpdir"));

        DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
        Date date = new Date();

        folder = new File(tmpDir, "BgenGenotypeDataTest_" + dateFormat.format(date));

        Runtime.getRuntime().addShutdownHook(new Thread() {
            @Override
            public void run() {
                System.out.println("Removing tmp dir and files");
                for (File file : folder.listFiles()) {
                    System.out.println(" - Deleting: " + file.getAbsolutePath());
                    file.delete();
                }
                System.out.println(" - Deleting: " + folder.getAbsolutePath());
                folder.delete();
            }
        });

        folder.mkdir();

        System.out.println("Temp folder with output of this test: " + folder.getAbsolutePath());
    }

//    @Test
//    public void testSomeMethod() {
//        // TODO review the generated test code and remove the default call to fail.
//    }

//    @Test
//    public void testBgenGenotypeData() throws URISyntaxException, IOException {
//        File bgenFile = getTestResourceFile("/bgenExamples/example.16bits.bgen");
//        Path target = Paths.get(folder.toString(), bgenFile.getName());
//        Files.copy(bgenFile.toPath(), target);
//        bgenGenotypeData = new BgenGenotypeData(target.toFile());
//        for (GeneticVariant variant : bgenGenotypeData) {
//            System.out.printf("%s %s %d %s%n", variant.getPrimaryVariantId(), variant.getSequenceName(), variant.getStartPos(), variant.getVariantAlleles());
//            variant.getSampleGenotypeProbilities();
//        }
//    }

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

            // Set the maximum error from the bit representation as shown within the BGEN specification.
            double maximumError = 1 / (Math.pow(2, bitRepresentation) - 1);
            // Replace the maximum error by 1*10^4 if the original maximum error was lower, as the precision of
            // the probabilities in the .gen file is often not greater than 4 decimals.
            maximumError = Math.max(maximumError, 1e-4);

            // Load the bgen file from a temporary folder
            Path bgenFile = Paths.get(folder.toString(), origBgenFile.getName());
            Files.copy(origBgenFile.toPath(), bgenFile);

            // The version 1.1 example file does not have sample identifiers
            if (origBgenFile.getName().equals("example.v11.bgen")) {
                bgenGenotypeData = new BgenGenotypeData(bgenFile.toFile(), exampleSampleFile);
                assertFalse(bgenGenotypeData.areSampleIdentifiersPresent());
                assertEquals(bgenGenotypeData.getSampleAnnotationsMap(), genGenotypeData.getSampleAnnotationsMap());
            } else {
                bgenGenotypeData = new BgenGenotypeData(bgenFile.toFile());
                assertTrue(bgenGenotypeData.areSampleIdentifiersPresent());
                // Test the equality of sample names and sequence names
                assertEquals(bgenGenotypeData.getSampleNames(), genGenotypeData.getSampleNames());
                assertEquals(bgenGenotypeData.getSampleAnnotationsMap(), new LinkedHashMap<>());
            }

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

    static void assertProbabilityEquality(float[][] actual, float[][] expected, double maximumError) {
        assertEquals(actual.length, expected.length);
        for (int i = 0; i < actual.length; i++) {
            assertEquals(actual[i].length, expected[i].length,
                    String.format("number of probabilities not equal for sample %d", i));
            for (int y = 0; y < actual[i].length; y++) {
                assertEquals(actual[i][y], expected[i][y],
                        maximumError, String.format("prob %d in sample %d", y, i));
            }
        }
    }

    static void assertProbabilityEquality(double[][] actual, double[][] expected) {
        assertEquals(actual.length, expected.length);
        for (int i = 0; i < actual.length; i++) {
            assertEquals(actual[i].length, expected[i].length,
                    String.format("number of probabilities not equal for sample %d", i));
            for (int y = 0; y < actual[i].length; y++) {
                assertEquals(actual[i][y], expected[i][y], String.format("prob %d in sample %d", y, i));
            }
        }
    }

    @Test
    public void complexBgenGenotypeDataTest() throws IOException {

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
        List<double[][]> expectedBgenProbabilities = Arrays.asList(
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

        List<double[][][]> expectedPhasedBgenProbabilities = Arrays.asList(
                null,
                new double[][][]{{{1, 0}}, {{1, 0}}, {{1, 0}}, {{0, 1}}},
                new double[][][]{{{1, 0}}, {{0, 1}, {0, 1}}, {{1, 0}, {1, 0}}, {{1, 0}, {0, 1}}},
                null,
                new double[][][]{{{1, 0}}, {{1, 0}, {1, 0}, {1, 0}}, {{1, 0}, {1, 0}, {0, 1}}, {{1, 0}, {0, 1}}},
                new double[][][]{{{1, 0, 0, 0}}, {{0, 1, 0, 0}}, {{0, 0, 1, 0}}, {{0, 0, 0, 1}}});

        // Loop through the complex bgen files and load the corresponding copy from a temporary folder
        for (File origBgenFile : complexFiles) {
            // Copy the BGEN file and load this copy
            Path bgenFile = Paths.get(folder.toString(), origBgenFile.getName());
            Files.copy(origBgenFile.toPath(), bgenFile);
            bgenGenotypeData = new BgenGenotypeData(bgenFile.toFile(), exampleComplexSampleFile);

            // Test the equality of sample names and sequence names
            assertEquals(bgenGenotypeData.getSampleNames(), expectedSampleNames);
            assertEquals(bgenGenotypeData.getSeqNames(), Collections.singletonList("01"));

            assertEquals(bgenGenotypeData.getVariantIdMap().keySet(), new HashSet<>(expectedVariantIds));
            // Loop through variants and check their similarity.
            int variantIndex = 0;
            for (GeneticVariant bgenVariant : bgenGenotypeData) {
                // Check if these variants are equal:
                assertEquals(bgenVariant.getPrimaryVariantId(), expectedVariantIds.get(variantIndex));
                assertEquals(bgenVariant.getStartPos(), expectedVariantPositions.get(variantIndex).intValue());
                assertEquals(bgenVariant.getAlternativeAlleles(), alleles.get(variantIndex));

                // V2 has an alternative id, the other variants don't. Check this.
                if (bgenVariant.getPrimaryVariantId().equals("V2")) {
                    assertEquals(bgenVariant.getAlternativeVariantIds().get(0), "V2.1");
                } else {
                    assertEquals(bgenVariant.getAlternativeVariantIds().size(), 0);
                }

                // Check the equality of probabilities.
                // First check if the bgenProbabilities are according to the expected stuff
                if (Arrays.asList(0, 1, 2, 3, 4, 8, 9).contains(variantIndex)) {
                    double[][] bgenProbabilities = bgenVariant.getSampleGenotypeProbabilitiesComplex();
                    assertProbabilityEquality(bgenProbabilities, expectedBgenProbabilities.get(variantIndex));
                }

                // Check if the regular probabilities are according to the expected stuff
                // (a lot should be coded as being missing)
                if (Arrays.asList(0, 3, 4).contains(variantIndex)) {
                    float[][] probabilities = bgenVariant.getSampleGenotypeProbilities();
                    assertProbabilityEquality(
                            probabilities,
                            expectedProbabilities.get(variantIndex),
                            0); // Maximum error is 0 as no decimal values are expected
                }

                // Check if the phased probabilities are correct as well
                if (Arrays.asList(1, 2, 4, 5).contains(variantIndex)) {
                    // First we have to check if the variant is phased
                    assertFalse(bgenVariant.getSamplePhasing().contains(false));
                    // Only then get the phased probabilities
                    double[][][] phasedBgenProbabilities = bgenVariant.getSampleGenotypeProbabilitiesPhased();
                    assertTrue(Arrays.deepEquals(phasedBgenProbabilities, expectedPhasedBgenProbabilities.get(variantIndex)));
                } else if (Arrays.asList(6, 7).contains(variantIndex)) {
                    // These are also phased, but probabilities are not interesting enough...
                    assertFalse(bgenVariant.getSamplePhasing().contains(false));
                } else {
                    // These are not phased, which we have to test for
                    assertFalse(bgenVariant.getSamplePhasing().contains(true));
                    // Calling the method then should generate an exception
                    try {
                        bgenVariant.getSampleGenotypeProbabilitiesPhased();
                        fail("bgenVariant.getSampleGenotypeProbabilitiesBgenPhased() did not raise a " +
                                "GenotypeDataException while phased data was not available");
                    } catch (GenotypeDataException e) {
                        assertEquals(e.getMessage(), "Phased data not available");
                    }
                }
                variantIndex++;
            }

        }
    }

    @Test
    public void testReaderHaplotypes() throws URISyntaxException, IOException {
        // Get the bgen input file to test with
        File bgenFile = getTestResourceFile("/bgenExamples/haplotypes.bgen");
        // Get the corresponding haplotype & sample files
        File hapsFile = getTestResourceFile("/bgenExamples/genFiles/haplotypes.haps");
        File oxfordSampleFile = getTestResourceFile("/bgenExamples/genFiles/haplotypes.sample");

        // Copy the input file
        Path target = Paths.get(folder.toString(), bgenFile.getName());
        Files.copy(bgenFile.toPath(), target);
        // Load the haps and bgen data
        HapsGenotypeData hapsGenotypeData = new HapsGenotypeData(hapsFile, oxfordSampleFile);
        bgenGenotypeData = new BgenGenotypeData(target.toFile(), exampleComplexSampleFile);

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
