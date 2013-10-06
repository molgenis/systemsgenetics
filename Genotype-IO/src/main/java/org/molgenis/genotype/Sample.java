package org.molgenis.genotype;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import org.molgenis.genotype.annotation.CaseControlAnnotation;
import org.molgenis.genotype.annotation.SexAnnotation;

public class Sample {

    private final String id;
    private final String familyId;
    private final Map<String, Object> annotationValues;

    public Sample(String id, String familyId, Map<String, Object> annotationValues) {
		
		if(id == null){
			throw new GenotypeDataException("Cannot create sample with a null ID");
		}
		
        this.id = id;
        this.familyId = familyId;
        if (annotationValues == null) {
            this.annotationValues = new HashMap<String, Object>();
        } else {
            this.annotationValues = annotationValues;
        }

    }

    public String getId() {
        return id;
    }

    public String getFamilyId() {
        return familyId;
    }

    public void putAnnotationValues(String annotationName, Object value) {
        annotationValues.put(annotationName, value);
    }

    public Map<String, ?> getAnnotationValues() {
        return Collections.unmodifiableMap(annotationValues);
    }

    public SexAnnotation getSex() {
        if (annotationValues.containsKey(GenotypeData.SEX_SAMPLE_ANNOTATION_NAME)) {
            return (SexAnnotation) annotationValues.get(GenotypeData.SEX_SAMPLE_ANNOTATION_NAME);
        } else {
            return SexAnnotation.UNKNOWN;
        }
    }
	
	public boolean isIncluded() {
		if(annotationValues.containsKey(GenotypeData.BOOL_INCLUDE_SAMPLE)){
			return (Boolean) annotationValues.get(GenotypeData.BOOL_INCLUDE_SAMPLE);
		} else {
			return true;
		}
	}
	
	public CaseControlAnnotation getCaseControlAnnotation() {
		if (annotationValues.containsKey(GenotypeData.CASE_CONTROL_SAMPLE_ANNOTATION_NAME)) {
            return (CaseControlAnnotation) annotationValues.get(GenotypeData.CASE_CONTROL_SAMPLE_ANNOTATION_NAME);
        } else {
            return CaseControlAnnotation.UNKNOWN;
        }
	}
	
	/**
	 * Get the ID of the father. "0" if not set. 
	 *
	 * @return 
	 */
	public String getFatherId(){
		if(annotationValues.containsKey(GenotypeData.FATHER_SAMPLE_ANNOTATION_NAME)){
			return (String) annotationValues.get(GenotypeData.FATHER_SAMPLE_ANNOTATION_NAME);
		} else {
			return "0";
		}
	}

	/**
	 * Get the ID of the father. "0" if not set. 
	 *
	 * @return 
	 */
	public String getMotherId(){
		if(annotationValues.containsKey(GenotypeData.MOTHER_SAMPLE_ANNOTATION_NAME)){
			return (String) annotationValues.get(GenotypeData.MOTHER_SAMPLE_ANNOTATION_NAME);
		} else {
			return "0";
		}
	}
	
	public float getMissingRate(){
		if(annotationValues.containsKey(GenotypeData.SAMPLE_MISSING_RATE_FLOAT)){
			return (Float) annotationValues.get(GenotypeData.SAMPLE_MISSING_RATE_FLOAT);
		} else {
			return Float.NaN;
		}
	}
	
    @Override
    public String toString() {
        return "Sample [id=" + id + ", familyId=" + familyId + ", annotationValues=" + annotationValues + "]";
    }

	@Override
	public int hashCode() {
		return this.id.hashCode();
	}

	@Override
	public boolean equals(Object obj) {
		if (obj == null) {
			return false;
		}
		if (getClass() != obj.getClass()) {
			return false;
		}
		final Sample other = (Sample) obj;
		if ((this.id == null) ? (other.id != null) : !this.id.equals(other.id)) {
			return false;
		}
		if ((this.familyId == null) ? (other.familyId != null) : !this.familyId.equals(other.familyId)) {
			return false;
		}
		return true;
	}

	
 
}
