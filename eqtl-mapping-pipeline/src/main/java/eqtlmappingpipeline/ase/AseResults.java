/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.ase;

import gnu.trove.map.hash.TIntObjectHashMap;
import java.util.HashMap;
import java.util.Iterator;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.variant.id.GeneticVariantId;
import umcg.genetica.collections.ChrPosMap;

/**
 *
 * @author Patrick Deelen
 */
public class AseResults implements Iterable<AseVariantAppendable> {

	private final ChrPosMap<AseVariantAppendable> results; //Empty chr hashmaps will crash the iterator
	private boolean encounteredBaseQuality = false;

	public AseResults() {
		results = new ChrPosMap<AseVariantAppendable>();
	}

	public synchronized void addResult(String chr, int pos, GeneticVariantId id, Allele a1, Allele a2, int a1Count, int a2Count, String sampleId, double a1MeanBaseQuality, double a2MeanBaseQuality) {
		
		addToResults(chr, pos, id, a1, a2, a1Count, a2Count, sampleId, a1MeanBaseQuality, a2MeanBaseQuality);
		encounteredBaseQuality = true;
	}
	
	public synchronized void addResult(String chr, int pos, GeneticVariantId id, Allele a1, Allele a2, int a1Count, int a2Count, String sampleId) {

		addToResults(chr, pos, id, a1, a2, a1Count, a2Count, sampleId, Double.NaN, Double.NaN);

	}
	
	private synchronized void addToResults(String chr, int pos, GeneticVariantId id, Allele a1, Allele a2, int a1Count, int a2Count, String sampleId, double a1MeanBaseQuality, double a2MeanBaseQuality) {

		AseVariantAppendable aseVariant = results.get(chr, pos);
		if (aseVariant == null) {
			aseVariant = new AseVariantAppendable(chr, pos, id, a1, a2);
			results.put(chr, pos, aseVariant);
		}
		aseVariant.addCounts(a1Count, a2Count, sampleId, a1MeanBaseQuality, a2MeanBaseQuality);
		
	}
	
	public Iterator<AseVariantAppendable> chrIterator(String chr){
		return results.chrIterator(chr);
	}

	@Override
	public Iterator<AseVariantAppendable> iterator() {
		return results.iterator();
	}
	
	public boolean isEncounteredBaseQuality() {
		return encounteredBaseQuality;
	}
	
	public int getCount(){
		return results.size();
	}
	
}