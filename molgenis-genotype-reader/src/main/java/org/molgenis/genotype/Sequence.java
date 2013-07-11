package org.molgenis.genotype;

import org.molgenis.genotype.variant.GeneticVariant;

/**
 * Represents a genetic sequence for example a chromosome
 * 
 * @author erwin
 * 
 */
public interface Sequence extends Iterable<GeneticVariant>
{
	String getName();

	Integer getLength();

	/**
	 * Returns if this sequence represents a chromosome or not
	 */
	boolean isChromosome();

	boolean isAutosome();

}
