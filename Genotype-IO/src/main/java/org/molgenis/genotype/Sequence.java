package org.molgenis.genotype;

/**
 * Represents a genetic sequence for example a chromosome
 * 
 * @author erwin
 * 
 */
public interface Sequence
{
	String getName();

	Integer getLength();

	/**
	 * Returns if this sequence represents a chromosome or not
	 */
	boolean isChromosome();

	boolean isAutosome();

}
