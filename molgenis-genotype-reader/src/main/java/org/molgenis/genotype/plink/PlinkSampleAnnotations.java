package org.molgenis.genotype.plink;

import java.util.LinkedHashMap;
import java.util.Map;

import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.annotation.SampleAnnotation.SampleAnnotationType;

public class PlinkSampleAnnotations
{

	public static Map<String, SampleAnnotation> getSampleAnnotations()
	{
		Map<String, SampleAnnotation> sampleAnnotations = new LinkedHashMap<String, SampleAnnotation>();

		SampleAnnotation fatherAnno = new SampleAnnotation(GenotypeData.FATHER_SAMPLE_ANNOTATION_NAME,
				GenotypeData.FATHER_SAMPLE_ANNOTATION_NAME, "", Annotation.Type.STRING, SampleAnnotationType.COVARIATE,
				false);
		sampleAnnotations.put(fatherAnno.getId(), fatherAnno);

		SampleAnnotation motherAnno = new SampleAnnotation(GenotypeData.MOTHER_SAMPLE_ANNOTATION_NAME,
				GenotypeData.MOTHER_SAMPLE_ANNOTATION_NAME, "", Annotation.Type.STRING, SampleAnnotationType.COVARIATE,
				false);
		sampleAnnotations.put(motherAnno.getId(), motherAnno);

		SampleAnnotation sexAnno = new SampleAnnotation(GenotypeData.SEX_SAMPLE_ANNOTATION_NAME,
				GenotypeData.SEX_SAMPLE_ANNOTATION_NAME, "", Annotation.Type.SEX, SampleAnnotationType.COVARIATE, false);
		sampleAnnotations.put(sexAnno.getId(), sexAnno);

		SampleAnnotation phenoAnno = new SampleAnnotation(GenotypeData.DOUBLE_PHENOTYPE_SAMPLE_ANNOTATION_NAME,
				GenotypeData.DOUBLE_PHENOTYPE_SAMPLE_ANNOTATION_NAME, "", Annotation.Type.FLOAT,
				SampleAnnotationType.PHENOTYPE, false);
		sampleAnnotations.put(phenoAnno.getId(), phenoAnno);

		return sampleAnnotations;
	}

}
