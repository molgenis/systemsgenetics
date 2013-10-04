package org.molgenis.genotype;

import java.io.Closeable;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 * Interface that represents genomic data, can be backed by different data types
 * 
 * @author erwin
 * 
 */
public interface GenotypeData extends Iterable<GeneticVariant>, Closeable
{

	public static final String FATHER_SAMPLE_ANNOTATION_NAME = "fatherId";
	public static final String MOTHER_SAMPLE_ANNOTATION_NAME = "motherId";
	public static final String SEX_SAMPLE_ANNOTATION_NAME = "sex_generic";
	public static final String DOUBLE_PHENOTYPE_SAMPLE_ANNOTATION_NAME = "phenotype";
    public static final String BOOL_INCLUDE_SAMPLE = "sampleInclude";
    public static final String CASE_CONTROL_SAMPLE_ANNOTATION_NAME = "caseControl";
	public static final String SAMPLE_MISSING_RATE_FLOAT = "sampleMissingRateDouble";

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
	
	/**
	 * Get variant annotations by id
	 */
	Map<String, Annotation> getVariantAnnotationsMap();
	
	/**
	 * Get variant annotations by id
	 */
	Map<String, SampleAnnotation> getSampleAnnotationsMap();
	
	/**
	 * 
	 *  
	 * @return 
	 */
	String[] getSampleNames();
	
	/**
	 * Return true if data contains probabilities or called genotypes (0,1,2)
	 * It is not 100% save to convert dosage that are not exactly 0, 1 or 2.
	 * Datasets that potentially contain these dosages must return false
	 */
	boolean isOnlyContaingSaveProbabilityGenotypes();

}
