/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.variantFilter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.molgenis.genotype.GenotypeDataException;
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

    public VariantCombinedFilter() {
        this.variantFilters = new ArrayList<VariantFilter>();
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

    @Override
    public boolean doesIdPassFilter(String id) {
        for (VariantFilter filter : variantFilters) {
            if (!filter.doesIdPassFilter(id)) {
                return false;
            }
        }
        return true;
    }

    public void add(VariantFilter v) {
        if (v == this) {
            throw new GenotypeDataException("Can not add a variant filter to it self");
        }
        variantFilters.add(v);
    }

}
