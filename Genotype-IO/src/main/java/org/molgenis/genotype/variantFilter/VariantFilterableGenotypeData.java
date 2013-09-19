/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.variantFilter;

/**
 *
 * @author Patrick Deelen
 */
public interface VariantFilterableGenotypeData {

	public void setFilter(VariantFilter filter);

	public VariantFilter getFilter();
}
