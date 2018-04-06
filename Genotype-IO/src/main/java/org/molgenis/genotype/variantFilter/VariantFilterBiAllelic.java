/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.variantFilter;

import org.molgenis.genotype.variant.GeneticVariant;

/**
 *
 * @author harmjan
 */
public class VariantFilterBiAllelic implements VariantFilter {

    @Override
    public boolean doesVariantPassFilter(GeneticVariant variant) {
        return variant.isBiallelic();
    }

    @Override
    public boolean doesIdPassFilter(String id) {
        return true;
    }

}
