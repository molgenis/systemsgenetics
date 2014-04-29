/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.ase;

import com.google.common.collect.Iterables;
import gnu.trove.map.hash.TIntObjectHashMap;
import java.util.HashMap;
import java.util.Iterator;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.variant.id.GeneticVariantId;

/**
 *
 * @author Patrick Deelen
 */
public class AseResults implements Iterable<AseVariant> {

	private final HashMap<String, TIntObjectHashMap<AseVariant>> results; //Empty chr hashmaps will crash the iterator

	public AseResults() {
		results = new HashMap<String, TIntObjectHashMap<AseVariant>>();
	}

	public synchronized void addResult(String chr, int pos, GeneticVariantId id, Allele a1, Allele a2, int a1Count, int a2Count) {

		TIntObjectHashMap<AseVariant> chrResults = results.get(chr);
		if (chrResults == null) {
			chrResults = new TIntObjectHashMap<AseVariant>();
			results.put(chr, chrResults);
		}

		AseVariant aseVariant = chrResults.get(pos);
		if (aseVariant == null) {
			aseVariant = new AseVariant(chr, pos, id, a1, a2);
			chrResults.put(pos, aseVariant);
		}
		aseVariant.addCounts(a1Count, a2Count);

	}

	@Override
	public Iterator<AseVariant> iterator() {
		
		return new AseResultIterator();
	}

	private class AseResultIterator implements Iterator<AseVariant> {

		Iterator<TIntObjectHashMap<AseVariant>> chrResultsIterator;
		Iterator<AseVariant> variantIterator;

		public AseResultIterator() {
			chrResultsIterator = results.values().iterator();

			if (chrResultsIterator.hasNext()) {
				variantIterator = chrResultsIterator.next().valueCollection().iterator();
			}
		}

		@Override
		public boolean hasNext() {
			if (variantIterator == null) {
				return chrResultsIterator.hasNext();
			}
			if (variantIterator.hasNext()) {
				return true;
			} else {
				return chrResultsIterator.hasNext();
			}
		}

		@Override
		public AseVariant next() {
			if (!variantIterator.hasNext()) {
				variantIterator = chrResultsIterator.next().valueCollection().iterator();
			}
			return variantIterator.next();

		}

		@Override
		public void remove() {
			variantIterator.remove();
		}
	}
	
	public int getCount(){
		int count = 0;
		
		for(TIntObjectHashMap<AseVariant> chrResults : results.values()){
			count += chrResults.size();
		}
		return count;
	}
	
}
