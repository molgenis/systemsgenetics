package org.molgenis.genotype;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.SampleAnnotation;

public abstract class AbstractGenotypeData implements GenotypeData
{
	@Override
	public List<Annotation> getVariantAnnotations()
	{
		return Collections.unmodifiableList(new ArrayList<Annotation>(getVariantAnnotationsMap().values()));
	}

	@Override
	public Annotation getVariantAnnotation(String annotationId)
	{
		return getVariantAnnotationsMap().get(annotationId);
	}

	@Override
	public List<SampleAnnotation> getSampleAnnotations()
	{
		return Collections.unmodifiableList(new ArrayList<SampleAnnotation>(getSampleAnnotationsMap().values()));
	}

	@Override
	public Annotation getSampleAnnotation(String annotationId)
	{
		return getVariantAnnotationsMap().get(annotationId);
	}

}
