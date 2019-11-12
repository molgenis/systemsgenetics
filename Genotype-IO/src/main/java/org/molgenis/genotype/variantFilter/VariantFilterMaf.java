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
public class VariantFilterMaf implements VariantFilter {

    private double maf;

    public VariantFilterMaf(double maf) {
        this.maf = maf;
	}

    public void setMafCutoff(float maf) {
        this.maf = maf;
    }

    @Override
    public boolean doesVariantPassFilter(GeneticVariant variant) {

		double observedMaf = variant.getMinorAlleleFrequency();
		
        if (Double.isNaN(observedMaf) || observedMaf < maf) {
            return false;
        }

        return true;

    }

    @Override
    public boolean doesIdPassFilter(String id) {
        return true;
    }
}
