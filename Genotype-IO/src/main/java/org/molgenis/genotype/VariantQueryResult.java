package org.molgenis.genotype;

import java.io.Closeable;

import org.molgenis.genotype.variant.GeneticVariant;

/**
 * The result of a VariantQuery. Must be closed after use
 * 
 * @author erwin
 * 
 */
public interface VariantQueryResult extends Closeable, Iterable<GeneticVariant>
{
}
