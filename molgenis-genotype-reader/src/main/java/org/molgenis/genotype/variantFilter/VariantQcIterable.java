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
public class VariantQcIterable implements Iterable<GeneticVariant> {

	private final Iterator<GeneticVariant> iterator;

	public VariantQcIterable(Iterator<GeneticVariant> iterator) {
		this.iterator = iterator;
	}
	
	@Override
	public Iterator<GeneticVariant> iterator() {
		return iterator;
	}
	
}
