package org.molgenis.genotype.variant;

import java.util.Collections;
import java.util.Map;
import org.molgenis.genotype.Alleles;

/**
 *
 * @author Patrick Deelen
 */
public class GenotypeRecordMap implements GenotypeRecord {

	private final Map<String, Object> fields;

	public GenotypeRecordMap(Map<String, Object> fields) {

		this.fields = Collections.unmodifiableMap(fields);

	}

	@Override
	public Object getGenotypeRecordData(String recordId) {

		return fields.get(recordId);

	}

	@Override
	public Alleles getSampleAlleles() {
		
		return (Alleles) fields.get("GT");
				
	}

	@Override
	public float[] getSampleProbs() {
		return (float[]) fields.get("GP");
	}

	@Override
	public float getSampleDosage() {
		return (Float) fields.get("DS");
	}

	@Override
	public boolean containsGenotypeRecord(String recordId) {
		return fields.containsKey(recordId);
	}
	
	
}
