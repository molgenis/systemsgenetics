package org.molgenis.genotype.multipart;

import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.ResourceTest;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.annotation.VcfAnnotation;
import org.molgenis.genotype.oxford.GenGenotypeData;
import org.molgenis.genotype.vcf.VcfGenotypeData;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.IOException;
import java.net.URISyntaxException;
import java.util.*;

import static org.testng.Assert.*;

public class MultiPartGenotypeDataTest extends ResourceTest {
    private RandomAccessGenotypeData chr1VarAnnotated;
    private RandomAccessGenotypeData chr2VarAnnotated;
    private RandomAccessGenotypeData chr1SampleAnnotated;
    private RandomAccessGenotypeData chr2SampleAnnotated;
    private RandomAccessGenotypeData chr2sampleAnnotatedWithWrongOrder;
    private RandomAccessGenotypeData multiPartVarAnnotated;
    private RandomAccessGenotypeData multiPartSampleAnnotated;
    private Map<String, Annotation> expectedVariantAnnotationMap;
    private Map<String, Annotation> expectedSampleAnnotationMap;

    @BeforeClass
    public void beforeClass() throws IOException, URISyntaxException
    {
        // Initialize maps of expected variant annotations and sample
        // annotations given compatible genotype data
        createExpectedVariantAnnotationMap();
        createExpectedSampleAnnotationMap();

        // Load the genotype data with variant annotations for two chromosomes
        chr1VarAnnotated = new VcfGenotypeData(
                getTestResourceFile("/multiPart/chr1.vcf.gz"),
                getTestResourceFile("/multiPart/chr1.vcf.gz.tbi"), 0.8);
        chr2VarAnnotated = new VcfGenotypeData(
                getTestResourceFile("/multiPart/chr2.vcf.gz"),
                getTestResourceFile("/multiPart/chr2.vcf.gz.tbi"), 0.8);

        // Load the genotype data with corresponding sample annotations
        chr1SampleAnnotated = new GenGenotypeData(
                getTestResourceFile("/multiPart/chr1.gen"),
                getTestResourceFile("/multiPart/all.sample"));
        chr2SampleAnnotated = new GenGenotypeData(
                getTestResourceFile("/multiPart/chr2.gen"),
                getTestResourceFile("/multiPart/all.sample"));

        // Load separate genotype data with sample annotations in a different order
        chr2sampleAnnotatedWithWrongOrder = new GenGenotypeData(
                getTestResourceFile("/multiPart/chr2.gen"),
                getTestResourceFile("/multiPart/all_unordered.sample"));

        // Create multipart genotype data
        multiPartVarAnnotated = new MultiPartGenotypeData(chr1VarAnnotated, chr2VarAnnotated);
        multiPartSampleAnnotated = new MultiPartGenotypeData(chr1SampleAnnotated, chr2SampleAnnotated);
    }

    /**
     * Creates a map of sample annotations.
     */
    private void createExpectedSampleAnnotationMap() {
        expectedSampleAnnotationMap = new LinkedHashMap<>();
        // Add phenotype annotations
        expectedSampleAnnotationMap.put("sampleMissingRateFloat",
                new SampleAnnotation("sampleMissingRateFloat", "missing",
                        "Missing data proportion of each individual", Annotation.Type.FLOAT,
                        SampleAnnotation.SampleAnnotationType.OTHER, false));
        expectedSampleAnnotationMap.put("pheno",
                new SampleAnnotation("pheno", "pheno", "", Annotation.Type.FLOAT,
                        SampleAnnotation.SampleAnnotationType.OTHER, false));
    }

    /**
     * Creates a map of variant annotations.
     */
    private void createExpectedVariantAnnotationMap() {
        expectedVariantAnnotationMap = new LinkedHashMap<>();
        expectedVariantAnnotationMap.put("PR",
                new VcfAnnotation("PR",
                        "Provisional reference allele, may not be based on real reference genome",
                        Annotation.Type.BOOLEAN, 0, false, false, false));
        expectedVariantAnnotationMap.put("AF",
                new VcfAnnotation("AF", "Allele frequency", Annotation.Type.FLOAT,
                null, false, true, false));
    }

    /**
     * Test that expects an incompatible multipart genotype data exception.
     */
    @Test(
            expectedExceptions = IncompatibleMultiPartGenotypeDataException.class,
            expectedExceptionsMessageRegExp =
                    "Incompatible multi part genotype data\\. " +
                            "All files should contain identical samples in same order\\. Found sample: Sample " +
                            "\\[id=s_1, familyId=s_1, annotationValues=\\{sampleMissingRateFloat=0\\.0, pheno=160\\.0\\}\\] " +
                            "expected: Sample " +
                            "\\[id=s_1, familyId=null, annotationValues=\\{\\}\\] for chr: 2"
    )
    public void testIncompatibleMultiPartGenotypeData() {
        new MultiPartGenotypeData(chr1VarAnnotated, chr2SampleAnnotated);
    }

    /**
     * Test that expects an incompatible multipart genotype data exception.
     */
    @Test(
            expectedExceptions = IncompatibleMultiPartGenotypeDataException.class,
            expectedExceptionsMessageRegExp =
                    "Incompatible multi part genotype data\\. " +
                            "All files should contain identical samples in same order\\. Found sample: Sample " +
                            "\\[id=s_2, familyId=s_2, annotationValues=\\{sampleMissingRateFloat=0\\.0, pheno=164\\.0\\}\\] " +
                            "expected: Sample " +
                            "\\[id=s_1, familyId=s_1, annotationValues=\\{sampleMissingRateFloat=0\\.0, pheno=160\\.0\\}\\] " +
                            "for chr: 2"
    )
    public void testIncompatibleMultiPartGenotypeDataUnordered() {
        new MultiPartGenotypeData(chr1SampleAnnotated, chr2sampleAnnotatedWithWrongOrder);
    }

    @Test
    public void testGetVariantAnnotations() {
        List<Annotation> expected = new ArrayList<>(expectedVariantAnnotationMap.values());
        assertEquals(chr1VarAnnotated.getVariantAnnotations(), expected);
        assertEquals(chr2VarAnnotated.getVariantAnnotations(), expected);
        assertEquals(multiPartVarAnnotated.getVariantAnnotations(), expected);
    }

    @Test
    public void testGetVariantAnnotation() {
        assertEquals(multiPartVarAnnotated.getVariantAnnotation("AF"),
                expectedVariantAnnotationMap.get("AF"));
        assertEquals(multiPartVarAnnotated.getVariantAnnotation("PR"),
                expectedVariantAnnotationMap.get("PR"));
    }

    @Test
    public void testGetVariantAnnotationsMap() {
        assertEquals(chr1VarAnnotated.getVariantAnnotationsMap(), expectedVariantAnnotationMap);
        assertEquals(chr2VarAnnotated.getVariantAnnotationsMap(), expectedVariantAnnotationMap);
        assertEquals(multiPartVarAnnotated.getVariantAnnotationsMap(), expectedVariantAnnotationMap);
    }

    @Test
    public void testGetSampleAnnotations() {
        List<Annotation> expected = new ArrayList<>(expectedSampleAnnotationMap.values());
        assertEquals(chr1SampleAnnotated.getSampleAnnotations(), expected);
        assertEquals(chr1SampleAnnotated.getSampleAnnotations(), expected);
        assertEquals(multiPartSampleAnnotated.getSampleAnnotations(), expected);
    }

    @Test
    public void testGetSampleAnnotation() {
        assertEquals(multiPartSampleAnnotated.getSampleAnnotation("pheno"),
                expectedSampleAnnotationMap.get("pheno"));
        assertEquals(multiPartSampleAnnotated.getSampleAnnotation("sampleMissingRateFloat"),
                expectedSampleAnnotationMap.get("sampleMissingRateFloat"));
     }

    @Test
    public void testGetSampleAnnotationsMap() {
        assertEquals(chr1SampleAnnotated.getSampleAnnotationsMap(), expectedSampleAnnotationMap);
        assertEquals(chr2SampleAnnotated.getSampleAnnotationsMap(), expectedSampleAnnotationMap);
        assertEquals(multiPartSampleAnnotated.getSampleAnnotationsMap(), expectedSampleAnnotationMap);
    }
}