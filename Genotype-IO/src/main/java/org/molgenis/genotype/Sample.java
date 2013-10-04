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

    /*
     * (non-Javadoc)
     * 
     * @see java.lang.Object#hashCode()
     */
    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + ((annotationValues == null) ? 0 : annotationValues.hashCode());
        result = prime * result + ((familyId == null) ? 0 : familyId.hashCode());
        result = prime * result + ((id == null) ? 0 : id.hashCode());
        return result;
    }

    /*
     * (non-Javadoc)
     * 
     * @see java.lang.Object#equals(java.lang.Object)
     */
    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        Sample other = (Sample) obj;
        if (annotationValues == null) {
            if (other.annotationValues != null) {
                return false;
            }
        } else if (!annotationValues.equals(other.annotationValues)) {
            return false;
        }
        if (familyId == null) {
            if (other.familyId != null) {
                return false;
            }
        } else if (!familyId.equals(other.familyId)) {
            return false;
        }
        if (id == null) {
            if (other.id != null) {
                return false;
            }
        } else if (!id.equals(other.id)) {
            return false;
        }
        return true;
    }
}
