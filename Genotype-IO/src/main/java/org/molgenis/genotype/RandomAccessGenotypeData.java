package org.molgenis.genotype;

import java.util.HashMap;
import java.util.List;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variantFilter.VariantFilter;


public interface RandomAccessGenotypeData extends GenotypeData {

	/**
	 * Get all sequence names
	 *
	 * @return List of String
	 */
	List<String> getSeqNames();

	/**
	 * Get all sequences in the data
	 *
	 * @return Iterable of Sequences
	 */
	Iterable<Sequence> getSequences();

	/**
	 * Get a Sequence buy it's name. Name is case sensitive
	 *
	 * @param name
	 * @return the Sequence or null if not found
	 */
	Sequence getSequenceByName(String name);

	/**
	 * Get the variants at the specified position
	 *
	 * @param seqName
	 * @param startPos
	 * @return all variants found at startPos, can be empty if none found
	 */
	Iterable<GeneticVariant> getVariantsByPos(String seqName, int startPos);

	/**
	 * Get the SNP variant at the specified position. Only one SNP possible per
	 * position.
	 *
	 * @param seqName
	 * @param startPos
	 * @return The SNP found at this startPos, will be null if not present
	 */
	GeneticVariant getSnpVariantByPos(String seqName, int startPos);

	/**
	 * Get all variants from a sequence
	 *
	 * @param seqName
	 * @return
	 */
	Iterable<GeneticVariant> getSequenceGeneticVariants(String seqName);

	/**
	 * Get all variants within the specified range
	 *
	 * @param seqName
	 * @param rangeStart start of range, inclusive
	 * @param rangeEnd end of range exclusive
	 * @return
	 */
	Iterable<GeneticVariant> getVariantsByRange(String seqName, int rangeStart, int rangeEnd);
	
	/**
	 * Get a HashMap with the variants that have a primairy ID.
	 * 
	 * Variants without an ID will always be ignored.
	 * 
	 * If multiple variants have the same ID an arbritary variant is selected.
	 * 
	 * It is possible to supply a variant filter to select a subset of variants
	 * 
	 * @param filter
	 * @return 
	 */
	HashMap<String, GeneticVariant> getVariantIdMap(VariantFilter filter);
	
	/**
	 * Variant ID map without filter is saved as cache, use this function to clear this cache.
	 */
	void clearVariantIdMap();
	
	/**
	 * Get a HashMap with the variants that have a primairy ID.
	 * 
	 * Variants without an ID will always be ignored.
	 * 
	 * If multiple variants have the same ID an arbritary variant is selected.
	 * 
	 * @param 
	 * @return 
	 */
	HashMap<String, GeneticVariant> getVariantIdMap();
}
