/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.variantFilter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 *
 * @author Patrick Deelen
 */
public class VariantCombinedFilter implements VariantFilter {

	private final List<VariantFilter> variantFilters;

	/**
	 * Allows specification of arbitrary number of variant filters
	 *
	 * Example: new VariantCombinedFilter(new VariantQcChecker(0.05, 0.95,
	 * 0.005), new VariantIdIncludeFilter(null));
	 *
	 * @param variantFilters
	 */
	public VariantCombinedFilter(VariantFilter... variantFilters) {
		this.variantFilters = Arrays.asList(variantFilters);
	}

	public VariantCombinedFilter(List<VariantFilter> variantFilters) {
		this.variantFilters = variantFilters == null ? new ArrayList<VariantFilter>() : variantFilters;
	}

	@Override
	public boolean doesVariantPassFilter(GeneticVariant variant) {

		for (VariantFilter filter : variantFilters) {
			if (!filter.doesVariantPassFilter(variant)) {
				return false;
			}
		}
		return true;


	}
}
