package org.molgenis.genotype;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import static org.molgenis.genotype.GenotypeData.SAMPLE_MISSING_RATE_FLOAT;

import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.variant.GeneticVariant;

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
		return getSampleAnnotationsMap().get(annotationId);
	}

	@Override
	public String[] getSampleNames() {
		
		String[] sampleNames = new String[getSamples().size()];
		
		int i = 0;
		for(Sample sample : getSamples()){
			sampleNames[i] = sample.getId();
			++i;
		}
		
		return sampleNames;
	}
	
	
@Override
	public void addSampleAnnotation(SampleAnnotation sampleAnnotation) {
		getSampleAnnotationsMap().put(sampleAnnotation.getName(), sampleAnnotation);
	}

	@Override
	public boolean containsSampleAnnotation(String annotationName) {
		return getSampleAnnotationsMap().containsKey(annotationName);
	}
	
	
	

}
