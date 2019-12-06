package org.molgenis.genotype.bgen;

import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.ResourceTest;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.plink.BedBimFamGenotypeData;
import org.molgenis.genotype.util.GenotypeDataCompareTool;
import org.molgenis.genotype.variant.GeneticVariant;
import org.testng.annotations.BeforeClass;
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

import static org.molgenis.genotype.bgen.BgenGenotypeDataTest.assertProbabilityEquality;
import static org.testng.Assert.*;

public class BgenGenotypeWriterTest extends ResourceTest {

    private File tmpOutputFolder;
    private String fileSep = System.getProperty("file.separator");
    private File complexFile = getTestResourceFile("/bgenExamples/complex.31bits.bgen");
    private File exampleFile = getTestResourceFile("/bgenExamples/example.16bits.bgen");

    public BgenGenotypeWriterTest() throws URISyntaxException {
    }

    @BeforeClass
    public void setUp() throws IOException, URISyntaxException {

        File tmpDir = new File(System.getProperty("java.io.tmpdir"));

        DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
        Date date = new Date();

        tmpOutputFolder = new File(tmpDir, "GenotypeBgenTest_" + dateFormat.format(date));


        Runtime.getRuntime().addShutdownHook(new Thread() {
            @Override
            public void run() {
                System.out.println("Removing tmp dir and files");
                for (File file : tmpOutputFolder.listFiles()) {
                    System.out.println(" - Deleting: " + file.getAbsolutePath());
                    //file.delete();
                }
                System.out.println(" - Deleting: " + tmpOutputFolder.getAbsolutePath());
                tmpOutputFolder.delete();
            }
        });

        tmpOutputFolder.mkdir();
    }

    @Test
    public void write() throws IOException, URISyntaxException {

        GenotypeData genotypeData = new BedBimFamGenotypeData(getTestBed6(), getTestBim6(), getTestFam6(), 2);
        BgenGenotypeWriter writer = new BgenGenotypeWriter(genotypeData);

        writer.write(tmpOutputFolder.getAbsolutePath() + fileSep + "test");

        GenotypeData genotypeDataWritten = RandomAccessGenotypeDataReaderFormats.BGEN.createGenotypeData(
                tmpOutputFolder.getAbsolutePath() + fileSep + "test", 0);

        assertTrue(GenotypeDataCompareTool.same(genotypeData, genotypeDataWritten));
    }

    @Test
    public void writeComplex() throws IOException, URISyntaxException {

        // Load the bgen file from a temporary folder
        Path copiedComplexFile = Paths.get(tmpOutputFolder.toString(), complexFile.getName());
        Files.copy(complexFile.toPath(), copiedComplexFile);

        BgenGenotypeData expectedGenotypeData = new BgenGenotypeData(copiedComplexFile.toFile());
        BgenGenotypeWriter writer = new BgenGenotypeWriter(expectedGenotypeData);

        writer.write(tmpOutputFolder.getAbsolutePath() + fileSep + "testcomplex");

        BgenGenotypeData bgenGenotypeData = new BgenGenotypeData(new File(
                tmpOutputFolder.getAbsolutePath() + fileSep + "testcomplex.bgen"));

        // Test the equality of sample names and sequence names
        assertEquals(bgenGenotypeData.getSampleNames(), expectedGenotypeData.getSampleNames());
        assertEquals(bgenGenotypeData.getSeqNames(), expectedGenotypeData.getSeqNames());

        assertEquals(bgenGenotypeData.getVariantIdMap(), bgenGenotypeData.getVariantIdMap());
        // Loop through variants and check their similarity.
        Iterator<GeneticVariant> actualIterator = bgenGenotypeData.iterator();

        int variantIndex = 0;
        for (GeneticVariant expectedVariant : expectedGenotypeData) {
            assertTrue(actualIterator.hasNext());
            GeneticVariant bgenVariant = actualIterator.next();
            // Check if these variants are equal:
            assertEquals(bgenVariant.getPrimaryVariantId(), expectedVariant.getPrimaryVariantId());
            assertEquals(bgenVariant.getStartPos(), expectedVariant.getStartPos());
            assertEquals(bgenVariant.getAlternativeAlleles(), expectedVariant.getAlternativeAlleles());

            assertEquals(bgenVariant.getAlternativeVariantIds(),
                    expectedVariant.getAlternativeVariantIds());

            double[][] bgenProbabilities = bgenVariant.getSampleGenotypeProbabilitiesComplex();
            assertProbabilityEquality(bgenProbabilities, expectedVariant.getSampleGenotypeProbabilitiesComplex());

            // Check if the regular probabilities are according to the expected stuff
            // (a lot should be coded as being missing)
            float[][] probabilities = bgenVariant.getSampleGenotypeProbilities();
            assertProbabilityEquality(
                    probabilities,
                    expectedVariant.getSampleGenotypeProbilities(),
                    0); // Maximum error is 0 as no decimal values are expected

            // Check if the phased probabilities are correct as well
            if (Arrays.asList(1, 2, 4, 5, 6, 7).contains(variantIndex)) {
                // First we have to check if the variant is phased
                assertFalse(bgenVariant.getSamplePhasing().contains(false));
                // Only then get the phased probabilities
                assertTrue(Arrays.deepEquals(
                        bgenVariant.getSampleGenotypeProbabilitiesPhased(),
                        expectedVariant.getSampleGenotypeProbabilitiesPhased()));
            } else {
                // These are not phased, which we have to test for
                assertFalse(bgenVariant.getSamplePhasing().contains(true));
            }
            variantIndex++;
        }
    }

    @Test
    public void largeBgenGenotypeWriterTest() throws IOException {

        BgenGenotypeData expectedGenotypeData = new BgenGenotypeData(exampleFile);

        BgenGenotypeWriter writer = new BgenGenotypeWriter(expectedGenotypeData);

        writer.write(tmpOutputFolder.getAbsolutePath() + fileSep + "testlarge");

        BgenGenotypeData bgenGenotypeData = new BgenGenotypeData(new File(
                tmpOutputFolder.getAbsolutePath() + fileSep + "testlarge.bgen"));

        double maximumError = 1 / (Math.pow(2, 16) - 1);

        assertTrue(bgenGenotypeData.areSampleIdentifiersPresent());
        // Test the equality of sample names and sequence names
        assertEquals(bgenGenotypeData.getSampleNames(), expectedGenotypeData.getSampleNames());
        assertEquals(bgenGenotypeData.getSampleAnnotationsMap(), expectedGenotypeData.getSampleAnnotationsMap());

        assertEquals(bgenGenotypeData.getSeqNames(), expectedGenotypeData.getSeqNames());

        Iterator<GeneticVariant> bgenIterator = bgenGenotypeData.iterator();
        // Loop through variants and check their similarity.
        for (GeneticVariant expectedVariant : expectedGenotypeData) {

            // Check if the bgen iterator also has a next variant.
            assertTrue(bgenIterator.hasNext(), "bgenIterator is emptied while Oxford genIterator is not.");
            GeneticVariant bgenVariant = bgenIterator.next();

            // Check if these are equal:
            assertEquals(bgenVariant.getPrimaryVariantId(), expectedVariant.getPrimaryVariantId());
            assertEquals(bgenVariant.getVariantAlleles(), expectedVariant.getVariantAlleles());
            assertEquals(bgenVariant.getAlternativeAlleles(), expectedVariant.getAlternativeAlleles());
            assertEquals(bgenVariant.getRefAllele(), expectedVariant.getRefAllele());
            assertEquals(bgenVariant.getSequenceName(), expectedVariant.getSequenceName());
            assertEquals(bgenVariant.getStartPos(), expectedVariant.getStartPos());
            assertEquals(bgenVariant.getAlleleCount(), expectedVariant.getAlleleCount());

            // Check if the probabilities are similar enough to the probabilities within the .gen file
            float[][] bgenProbabilities = bgenVariant.getSampleGenotypeProbilities();

            float[][] expectedProbabilities = expectedVariant.getSampleGenotypeProbilities();
            assertProbabilityEquality(bgenProbabilities, expectedProbabilities, maximumError);
        }
    }
}