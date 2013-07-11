package org.molgenis.genotype.plink;

import java.util.LinkedHashMap;
import java.util.Map;

import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.annotation.SampleAnnotation.SampleAnnotationType;

public class PlinkSampleAnnotations
{

	public static Map<String, SampleAnnotation> getSampleAnnotations()
	{
		Map<String, SampleAnnotation> sampleAnnotations = new LinkedHashMap<String, SampleAnnotation>();

		SampleAnnotation fatherAnno = new SampleAnnotation(PedMapGenotypeData.FATHER_SAMPLE_ANNOTATION_NAME,
				PedMapGenotypeData.FATHER_SAMPLE_ANNOTATION_NAME, "", Annotation.Type.STRING,
				SampleAnnotationType.OTHER, false);
		sampleAnnotations.put(fatherAnno.getId(), fatherAnno);

		SampleAnnotation motherAnno = new SampleAnnotation(PedMapGenotypeData.MOTHER_SAMPLE_ANNOTATION_NAME,
				PedMapGenotypeData.MOTHER_SAMPLE_ANNOTATION_NAME, "", Annotation.Type.STRING,
				SampleAnnotationType.OTHER, false);
		sampleAnnotations.put(motherAnno.getId(), motherAnno);

		SampleAnnotation sexAnno = new SampleAnnotation(PedMapGenotypeData.SEX_SAMPLE_ANNOTATION_NAME,
				PedMapGenotypeData.SEX_SAMPLE_ANNOTATION_NAME, "", Annotation.Type.INTEGER, SampleAnnotationType.OTHER,
				false);
		sampleAnnotations.put(sexAnno.getId(), sexAnno);

		SampleAnnotation phenoAnno = new SampleAnnotation(PedMapGenotypeData.PHENOTYPE_SAMPLE_ANNOTATION_NAME,
				PedMapGenotypeData.PHENOTYPE_SAMPLE_ANNOTATION_NAME, "", Annotation.Type.FLOAT,
				SampleAnnotationType.PHENOTYPE, false);
		sampleAnnotations.put(phenoAnno.getId(), phenoAnno);

		return sampleAnnotations;
	}

}
