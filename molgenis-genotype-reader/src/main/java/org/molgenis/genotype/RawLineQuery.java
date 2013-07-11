package org.molgenis.genotype;

/**
 * Query that returns the raw lines of an inpute file
 * 
 * @author erwin
 * 
 */
public interface RawLineQuery
{
	RawLineQueryResult executeQuery(String sequence, int startPos);
}
