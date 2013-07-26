/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.variantFilter;

import java.util.ArrayList;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 *
 * @author Patrick Deelen
 */
public class VariantCombinedFilter implements VariantFilter{
	
	private final ArrayList<VariantFilter> variantFilters;

	public VariantCombinedFilter(ArrayList<VariantFilter> variantFilters) {
		this.variantFilters = variantFilters == null ? new ArrayList<VariantFilter>() : variantFilters;
	}

	@Override
	public boolean doesVariantPassFilter(GeneticVariant variant) {
		for(VariantFilter filter : variantFilters){
			if(! filter.doesVariantPassFilter(variant)){
				return false;
			}
		}
		return true;
	}
	
}
