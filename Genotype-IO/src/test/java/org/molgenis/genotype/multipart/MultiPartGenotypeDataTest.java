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
        createExpectedVariantAnnotationMap();
        createExpectedSampleAnnotationMap();

        chr1VarAnnotated = new VcfGenotypeData(
                getTestResourceFile("/multiPart/chr1.vcf.gz"),
                getTestResourceFile("/multiPart/chr1.vcf.gz.tbi"), 0.8);
        chr2VarAnnotated = new VcfGenotypeData(
                getTestResourceFile("/multiPart/chr2.vcf.gz"),
                getTestResourceFile("/multiPart/chr2.vcf.gz.tbi"), 0.8);

        chr1SampleAnnotated = new GenGenotypeData(
                getTestResourceFile("/multiPart/chr1.gen"),
                getTestResourceFile("/multiPart/all.sample"));
        chr2SampleAnnotated = new GenGenotypeData(
                getTestResourceFile("/multiPart/chr2.gen"),
                getTestResourceFile("/multiPart/all.sample"));

        chr2sampleAnnotatedWithWrongOrder = new GenGenotypeData(
                getTestResourceFile("/multiPart/chr2.gen"),
                getTestResourceFile("/multiPart/all_unordered.sample"));

        multiPartVarAnnotated = new MultiPartGenotypeData(chr1VarAnnotated, chr2VarAnnotated);
        multiPartSampleAnnotated = new MultiPartGenotypeData(chr1SampleAnnotated, chr2SampleAnnotated);
    }

    private void createExpectedSampleAnnotationMap() {
        expectedSampleAnnotationMap = new LinkedHashMap<>();
        expectedSampleAnnotationMap.put("pheno",
                new SampleAnnotation("pheno", "pheno", "", Annotation.Type.FLOAT,
                        SampleAnnotation.SampleAnnotationType.OTHER, false));
        expectedSampleAnnotationMap.put("sampleMissingRateFloat",
                new SampleAnnotation("sampleMissingRateFloat", "missing",
                        "Missing data proportion of each individual", Annotation.Type.FLOAT,
                        SampleAnnotation.SampleAnnotationType.OTHER, false));
    }

    private void createExpectedVariantAnnotationMap() {
        expectedVariantAnnotationMap = new LinkedHashMap<>();
        expectedVariantAnnotationMap.put("AF",
                new VcfAnnotation("AF", "Allele frequency", Annotation.Type.FLOAT,
                null, false, true, false));
        expectedVariantAnnotationMap.put("PR",
                new VcfAnnotation("PR",
                "Provisional reference allele, may not be based on real reference genome",
                Annotation.Type.BOOLEAN, 0, false, false, false));
    }

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
        assertTrue(chr1VarAnnotated.getVariantAnnotations().size() == expected.size() &&
                chr1VarAnnotated.getVariantAnnotations().containsAll(expected) &&
                expected.containsAll(chr1VarAnnotated.getVariantAnnotations()));
        assertTrue(chr2VarAnnotated.getVariantAnnotations().size() == expected.size() &&
                chr2VarAnnotated.getVariantAnnotations().containsAll(expected) &&
                expected.containsAll(chr2VarAnnotated.getVariantAnnotations()));
        assertTrue(multiPartVarAnnotated.getVariantAnnotations().size() == expected.size() &&
                multiPartVarAnnotated.getVariantAnnotations().containsAll(expected) &&
                expected.containsAll(multiPartVarAnnotated.getVariantAnnotations()));
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
        assertTrue(multiPartSampleAnnotated.getSampleAnnotations().size() == expected.size() &&
                multiPartSampleAnnotated.getSampleAnnotations().containsAll(expected) &&
                expected.containsAll(multiPartSampleAnnotated.getSampleAnnotations()));
    }

    @Test
    public void testGetSampleAnnotation() {
        assertEquals(multiPartVarAnnotated.getSampleAnnotation("pheno"),
                expectedSampleAnnotationMap.get("pheno"));
        assertEquals(multiPartVarAnnotated.getSampleAnnotation("sampleMissingRateFloat"),
                expectedSampleAnnotationMap.get("sampleMissingRateFloat"));
     }

    @Test
    public void testGetSampleAnnotationsMap() {
        assertEquals(chr1SampleAnnotated.getSampleAnnotationsMap(), expectedSampleAnnotationMap);
        assertEquals(chr2SampleAnnotated.getSampleAnnotationsMap(), expectedSampleAnnotationMap);
        assertEquals(multiPartSampleAnnotated.getSampleAnnotationsMap(), expectedSampleAnnotationMap);
    }
}