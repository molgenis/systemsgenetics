/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.util;

import com.google.common.base.Function;
import com.google.common.collect.Iterators;
import java.util.Iterator;
import java.util.List;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.variant.GenotypeRecord;
import org.molgenis.genotype.variant.GenotypeRecordGt;
import org.molgenis.genotype.variant.GenotypeRecordProb;

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
	
	public static FixedSizeIterable<GenotypeRecord> createIteratorFromProbs(float[][] probs){
		return new FixedSizeIterableProbs(probs);
	}
	
	private static class ProbsToProbsRecordFunction implements Function<float[], GenotypeRecord>{

		@Override
		public GenotypeRecord apply(float[] input) {
			return new GenotypeRecordProb(input);
		}

	}
	
	private static class FixedSizeIterableProbs implements FixedSizeIterable<GenotypeRecord>{

		private final float[][] probs;

		public FixedSizeIterableProbs(float[][] probs) {
			this.probs = probs;
		}
		
		@Override
		public int size() {
			return probs.length;
		}

		@Override
		public Iterator<GenotypeRecord> iterator() {
			return Iterators.transform(Iterators.forArray(probs), new ProbsToProbsRecordFunction());
		}
		
	}
	
}
