/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.variant;

import java.util.HashMap;
import org.molgenis.genotype.GenotypeDataException;

/**
 *
 * @author Patrick Deelen
 */
public class GenotypeRecordType {

	public enum Type {
		INTEGER, FLOAT, STRING, CHAR, ALLELES, INTEGER_ARRAY, FLOAT_ARRAY, STRING_ARRAY, CHAR_ARRAY;
	}
	private final String id;
	private final int number;
	private final Type type;
	private final String description;
	
	private static final HashMap<String, Type> RESEVERED_IDS;
	
	static{
		RESEVERED_IDS = new HashMap<String, Type>();
		RESEVERED_IDS.put("GT", Type.ALLELES);
	}

	public GenotypeRecordType(String id, int number, Type type, String description) {
		
		if(checkNotOverwriteReserved(id, type)){
			throw new GenotypeDataException("Expected for genotype record: " + id + " type " + RESEVERED_IDS.get(id).name() + " but found: " + type);
		}
		
		this.id = id;
		this.number = number;
		this.type = type;
		this.description = description;
	}

	public GenotypeRecordType(String id, int number, String typeName, String description) {
		this.id = id;
		this.number = number;
		try {
			this.type = Type.valueOf(typeName.toUpperCase());
		} catch (IllegalArgumentException ex) {
			throw new GenotypeDataException("Invalid genotype record type", ex);			
		}
		
		if(checkNotOverwriteReserved(id, type)){
			throw new GenotypeDataException("Expected for genotype record: " + id + " type " + RESEVERED_IDS.get(id).name() + " but found: " + type);
		}
		
		this.description = description;
	}
	
	private boolean checkNotOverwriteReserved(String id, Type type){
		if(RESEVERED_IDS.containsKey(id)){
			if(RESEVERED_IDS.get(id).equals(type)){
				return false;
			} else {
				return true;
			}
		} else {
			return false;
		}
	}
}
