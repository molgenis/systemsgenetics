/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.variantFilter;

import java.util.Iterator;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 *
 * @author Patrick Deelen
 */
public class VariantFilterIterable implements Iterable<GeneticVariant> {

	private final Iterator<GeneticVariant> iterator;

	public VariantFilterIterable(Iterator<GeneticVariant> iterator) {
		this.iterator = iterator;
	}
	
	@Override
	public Iterator<GeneticVariant> iterator() {
		return iterator;
	}
	
}
