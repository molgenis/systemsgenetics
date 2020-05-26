package org.molgenis.genotype.variantFilter;

import org.molgenis.genotype.variant.GeneticVariant;

/**
 *
 * @author Patrick Deelen
 */
public class VariantFilterMachR2 implements VariantFilter{

	private final double minMachR2;

	public VariantFilterMachR2(double minMachR2) {
		this.minMachR2 = minMachR2;
	}
	
	@Override
	public boolean doesVariantPassFilter(GeneticVariant variant) {
		return variant.getMachR2() >= minMachR2;
	}

	@Override
	public boolean doesIdPassFilter(String id) {
		return true;
	}

}
