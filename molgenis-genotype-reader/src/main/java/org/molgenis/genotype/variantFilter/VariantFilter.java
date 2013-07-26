/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.variantFilter;

import org.molgenis.genotype.variant.GeneticVariant;

/**
 *
 * @author Patrick Deelen
 */
public interface VariantFilter {
	
	public boolean doesVariantPassFilter(GeneticVariant variant);
	
	
}
