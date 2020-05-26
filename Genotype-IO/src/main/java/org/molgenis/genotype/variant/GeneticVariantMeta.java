package org.molgenis.genotype.variant;

public interface GeneticVariantMeta
{
	enum Type {INTEGER, FLOAT, STRING, CHAR, ALLELES, INTEGER_LIST, FLOAT_LIST, STRING_LIST, CHAR_LIST; }
	
	Iterable<String> getRecordIds();
	
	/**
	 * 
	 * @param recordId
	 * @return record data type or null if a record with the given id does not exist 
	 */
	Type getRecordType(String recordId);
}
