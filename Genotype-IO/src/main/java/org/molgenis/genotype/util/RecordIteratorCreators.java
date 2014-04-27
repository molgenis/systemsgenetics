/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.util;

import com.google.common.base.Function;
import com.google.common.collect.Iterators;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.variant.GenotypeRecord;
import org.molgenis.genotype.variant.GenotypeRecordGt;

/**
 *
 * @author Patrick Deelen
 */
public class RecordIteratorCreators {
	
	public static FixedSizeIterable<GenotypeRecord> createIteratorFromAlleles(List<Alleles> alleles){
		return new FixedSizeIterableAlleles(alleles);
	}
	
	private static class AllelesToAllelesRecordFunction implements Function<Alleles, GenotypeRecord>{

		@Override
		public GenotypeRecord apply(Alleles input) {
			return new GenotypeRecordGt(input);
		}

	}
	
	private static class FixedSizeIterableAlleles implements FixedSizeIterable<GenotypeRecord>{

		private final List<Alleles> alleles;

		public FixedSizeIterableAlleles(List<Alleles> alleles) {
			this.alleles = alleles;
		}
		
		@Override
		public int size() {
			return alleles.size();
		}

		@Override
		public Iterator<GenotypeRecord> iterator() {
			return Iterators.transform(alleles.iterator(), new AllelesToAllelesRecordFunction());
		}
		
	}
	
	
	
}
