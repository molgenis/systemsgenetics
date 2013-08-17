/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.variantFilter;

import java.util.HashSet;
import java.util.Set;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 *
 * @author Patrick Deelen
 */
public class VariantIdIncludeFilter implements VariantFilter{

	private final Set<String> include;
	
	public VariantIdIncludeFilter(Set<String> include) {
		this.include = include == null ? new HashSet<String>() : include;
	}

	@Override
	public boolean doesVariantPassFilter(GeneticVariant variant) {
		return include.contains(variant.getPrimaryVariantId());
	}
	
	public void addIdToInclude(String id){
		include.add(id);
	}
	
}
