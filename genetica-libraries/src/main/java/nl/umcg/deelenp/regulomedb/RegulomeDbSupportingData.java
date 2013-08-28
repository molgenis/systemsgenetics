/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.deelenp.regulomedb;

import java.util.regex.Pattern;

/**
 *
 * @author Patrick Deelen
 */
public class RegulomeDbSupportingData {
	
	private final String supportClass;
	private final String supportMethod;
	private final String supportValue;
	
	private static final Pattern PIPE_PATTERN = Pattern.compile("\\|");

	public RegulomeDbSupportingData(String supportData) throws Exception {
		String[] supportDataElements = PIPE_PATTERN.split(supportData);
		
		if(supportDataElements.length == 3){
			this.supportClass = supportDataElements[0];
			this.supportMethod = supportDataElements[1];
			this.supportValue = supportDataElements[2];
		} else if (supportDataElements.length == 2){
			this.supportClass = supportDataElements[0];
			this.supportMethod = supportDataElements[1];
			this.supportValue = "";
		} else if (supportDataElements.length == 4){
			this.supportClass = supportDataElements[0];
			this.supportMethod = supportDataElements[1];
			this.supportValue = supportDataElements[2] + "|" + supportDataElements[3];
		}else {
			throw new Exception("Error in RegulomeDB support data. Expected three elements separated by | found: " + supportData);
		}
				
		
	}

	public String getSupportClass() {
		return supportClass;
	}

	public String getSupportMethod() {
		return supportMethod;
	}

	public String getSupportValue() {
		return supportValue;
	}

	@Override
	public int hashCode() {
		int hash = 7;
		hash = 29 * hash + (this.supportClass != null ? this.supportClass.hashCode() : 0);
		hash = 29 * hash + (this.supportMethod != null ? this.supportMethod.hashCode() : 0);
		hash = 29 * hash + (this.supportValue != null ? this.supportValue.hashCode() : 0);
		return hash;
	}

	@Override
	public boolean equals(Object obj) {
		if (obj == null) {
			return false;
		}
		if (getClass() != obj.getClass()) {
			return false;
		}
		final RegulomeDbSupportingData other = (RegulomeDbSupportingData) obj;
		if ((this.supportClass == null) ? (other.supportClass != null) : !this.supportClass.equals(other.supportClass)) {
			return false;
		}
		if ((this.supportMethod == null) ? (other.supportMethod != null) : !this.supportMethod.equals(other.supportMethod)) {
			return false;
		}
		if ((this.supportValue == null) ? (other.supportValue != null) : !this.supportValue.equals(other.supportValue)) {
			return false;
		}
		return true;
	}
	
	
	
}
