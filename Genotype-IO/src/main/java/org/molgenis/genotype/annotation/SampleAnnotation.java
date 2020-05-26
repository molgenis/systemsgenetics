package org.molgenis.genotype.annotation;

public class SampleAnnotation extends Annotation
{
	public enum SampleAnnotationType
	{
		COVARIATE, PHENOTYPE, OTHER
	}

	private SampleAnnotationType sampleAnnotationType;
	private boolean list;

	public SampleAnnotation(String id, String name, String description, Type type,
			SampleAnnotationType sampleAnnotationType, boolean list)
	{
		super(id, name, description, type);
		this.sampleAnnotationType = sampleAnnotationType;
		this.list = list;
	}

	@Override
	public boolean isList()
	{
		return list;
	}

	public SampleAnnotationType getSampleAnnotationType()
	{
		return sampleAnnotationType;
	}

}
