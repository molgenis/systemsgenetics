/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.methylation;

import java.util.LinkedHashMap;

/**
 *
 * @author Marc Jan
 */
public class SoftfileAnnotation {

	private String title = null;
	private String accession = null;
	private String meshTerms = null;
	private LinkedHashMap<String, String> annotationInformation = new LinkedHashMap<String, String>();

	public SoftfileAnnotation() {
	}

	public String getAccession() {
		return accession;
	}

	public void setAccession(String accession) {
		this.accession = accession;
	}

	public void putAnnotationInformation(String key, String value) {
		this.annotationInformation.put(key, value);
	}

	public LinkedHashMap<String, String> getAnnotationInformation() {
		return annotationInformation;
	}

	public void setAnnotationInformation(LinkedHashMap<String, String> annotationInformation) {
		this.annotationInformation = annotationInformation;
	}

	public String getTitle() {
		return title;
	}

	public void setTitle(String title) {
		this.title = title;
	}

	public String getMeshTerms() {
		return meshTerms;
	}

	public void setMeshTerms(String meshTerms) {
		this.meshTerms = meshTerms;
	}
}
