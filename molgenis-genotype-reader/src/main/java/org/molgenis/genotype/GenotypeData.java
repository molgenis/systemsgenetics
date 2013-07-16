package org.molgenis.genotype;

import java.io.IOException;
import java.util.List;

import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 * Interface that represents genomic data, can be backed by different data types
 * 
 * @author erwin
 * 
 */
public interface GenotypeData extends Iterable<GeneticVariant>
{

	public static final String FATHER_SAMPLE_ANNOTATION_NAME = "fatherId";
	public static final String MOTHER_SAMPLE_ANNOTATION_NAME = "motherId";
	public static final String SEX_SAMPLE_ANNOTATION_NAME = "sex_generic";
	public static final String DOUBLE_PHENOTYPE_SAMPLE_ANNOTATION_NAME = "phenotype";

	/**
	 * Get all possible variant annotations
	 * 
	 * @return List of Annotation
	 */
	List<Annotation> getVariantAnnotations();

	/**
	 * Get a specific variant annotation
	 * 
	 * @param annotationId
	 * @return The Annotation or null if not found
	 * @throws IOException
	 */
	Annotation getVariantAnnotation(String annotationId);

	/**
	 * Get all posible sample annotations
	 * 
	 * @return
	 */
	List<SampleAnnotation> getSampleAnnotations();

	/**
	 * Get a specific sample annotation
	 * 
	 * @param annotationId
	 * @return The Annotation or null if not found
	 */
	Annotation getSampleAnnotation(String annotationId);

	/**
	 * 
	 * Get all samples
	 * 
	 * @return
	 */
	List<Sample> getSamples();

}
