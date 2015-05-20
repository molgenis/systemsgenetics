package org.molgenis.genotype;

import java.util.HashMap;
import java.util.Iterator;
import java.util.NoSuchElementException;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variantFilter.VariantFilter;

public abstract class AbstractRandomAccessGenotypeData extends AbstractGenotypeData implements RandomAccessGenotypeData {

	private HashMap<String, GeneticVariant> fullVariantMap = null;

	@Override
	public Sequence getSequenceByName(String name) {
		for (Sequence sequence : getSequences()) {
			if (sequence.getName().equals(name)) {
				return sequence;
			}
		}

		return null;
	}

	@Override
	public GeneticVariant getSnpVariantByPos(String seqName, int startPos) {
		Iterable<GeneticVariant> variants = getVariantsByPos(seqName, startPos);

		for (GeneticVariant variant : variants) {
			if (variant.isSnp()) {
				// only one SNP possible per position. Returning this SNP only
				return variant;
			}
		}

		return null;

	}

	@Override
	public HashMap<String, GeneticVariant> getVariantIdMap() {

		if (fullVariantMap == null) {
			fullVariantMap = getVariantIdMap(null);
		}
		return fullVariantMap;

	}
	
	@Override
	public void clearVariantIdMap(){
		fullVariantMap = null;
	}

	@Override
	public HashMap<String, GeneticVariant> getVariantIdMap(VariantFilter filter) {

		HashMap<String, GeneticVariant> variantIdMap = new HashMap<String, GeneticVariant>();

		for (GeneticVariant variant : this) {
			if (variant.getVariantId().getPrimairyId() != null && !variant.getPrimaryVariantId().equals("") && (filter == null || filter.doesVariantPassFilter(variant))) {
				variantIdMap.put(variant.getPrimaryVariantId(), variant);
			}
		}

		return variantIdMap;

	}

	@Override
	public Iterator<GeneticVariant> iterator() {
		return new GeneticVariantsIterator(this);
	}

	private static class GeneticVariantsIterator implements Iterator<GeneticVariant> {

		private Iterator<String> seqNames;
		private Iterator<GeneticVariant> seqGeneticVariants;
		private RandomAccessGenotypeData randomAccessGenotypeData;

		public GeneticVariantsIterator(RandomAccessGenotypeData randomAccessGenotypeData) {
			seqNames = randomAccessGenotypeData.getSeqNames().iterator();
			seqGeneticVariants = randomAccessGenotypeData.getSequenceGeneticVariants(seqNames.next()).iterator();
			this.randomAccessGenotypeData = randomAccessGenotypeData;
		}

		@Override
		public boolean hasNext() {
			return seqGeneticVariants.hasNext() || seqNames.hasNext();
		}

		@Override
		public GeneticVariant next() {
			if (seqGeneticVariants.hasNext()) {
				return seqGeneticVariants.next();
			}

			if (seqNames.hasNext()) {
				seqGeneticVariants = randomAccessGenotypeData.getSequenceGeneticVariants(seqNames.next()).iterator();
				return seqGeneticVariants.next();
			}

			throw new NoSuchElementException();
		}

		@Override
		public void remove() {
			throw new UnsupportedOperationException();
		}
	}
}
