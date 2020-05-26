package org.molgenis.genotype.variant.range;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.plink.BedBimFamGenotypeWriter;
import org.molgenis.genotype.util.ChromosomeComparator;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.GeneticVariantChrPosComparator;

/**
 * @author Patrick Deelen
 */
public class GeneticVariantRange implements Iterable<GeneticVariant> {

	private static final org.apache.log4j.Logger LOGGER = org.apache.log4j.Logger.getLogger(GeneticVariantRange.class);

	public static GeneticVariantRangeCreate createRangeFactory() {
		return new GeneticVariantRangeCreate();
	}

	public static GeneticVariantRangeCreate createRangeFactory(int expectedNumberVariants) {
		return new GeneticVariantRangeCreate(expectedNumberVariants);
	}

	private final List<GeneticVariant> variants;

	/**
	 * MAKE ABSOLUTELY SURE THAT VARIANTS IS SORTED BY CHR AND POS
	 *
	 * @param variants
	 */
	private GeneticVariantRange(List<GeneticVariant> variants) {
		this.variants = variants;
	}

	/**
	 * Get all variant from a specific seq/chr pos to (inclusive) a specified
	 * position
	 *
	 * @param startSeqName
	 * @param startPos
	 * @param stopSeqName
	 * @param stopPos
	 * @return
	 */
	public Iterable<GeneticVariant> getVariantsByRange(String startSeqName, int startPos, String stopSeqName, int stopPos) {

		if (ChromosomeComparator.chrASmallerChrB(startSeqName, stopSeqName)) {
			throw new GenotypeDataException("Cannot query range. Start chr is larger than stop chr");
		}

		if (startSeqName.equals(startSeqName) && startPos > stopPos) {
			throw new GenotypeDataException("Cannot query range. Start position is larger than stop postion and chr is identical");
		}

		//we need to find the first entry that is larger then startSeq and startPos.
		//aka the elements that matched the start search but whose left neighbour does not.
		//startIndex will be number of variants.size if no variants are within the range. The iteratble will handle this.
		int startIndex = binarySearchStartIndex(startSeqName, startPos, stopSeqName, stopPos);

		return new GeneticVariantRangeSubIterable(startIndex, stopSeqName, stopPos, variants);

	}

	/**
	 * WARNING: variants.size, end is handled by the iterable
	 *
	 * @param startSeqName
	 * @param startPos
	 * @param stopSeqName
	 * @param stopPos
	 * @return
	 */
	private int binarySearchStartIndex(String startSeqName, int startPos, String stopSeqName, int stopPos) {

		//define start search range
		int startSearch = 0;
		int endSearch = variants.size() - 1;

		while (endSearch >= startSearch) {

			int midSearch = startSearch + ((endSearch - startSearch) / 2);
			GeneticVariant midSearchVariant = variants.get(midSearch);

			//first test if current variants is larger than start of range. Then test if variant before is not in range
			if (ChromosomeComparator.chrASmallerChrB(midSearchVariant.getSequenceName(), startSeqName) || (midSearchVariant.getSequenceName().equals(startSeqName) && midSearchVariant.getStartPos() < startPos)) {
				//variant is smaller than start we have to search to right
				startSearch = midSearch + 1;
			} else if (ChromosomeComparator.chrASmallerEqualChrB(midSearchVariant.getSequenceName(), stopSeqName) && midSearchVariant.getStartPos() <= stopPos) {
				//variant is within the range of interest, lets check variant to the left

				if (midSearch == 0) {
					//We have found the start of our range. It streches from the beginning
					return 0;
				} else {

					GeneticVariant leftOfMidSearchVariant = variants.get(midSearch - 1);


					if (ChromosomeComparator.chrASmallerChrB(leftOfMidSearchVariant.getSequenceName(), startSeqName) || (leftOfMidSearchVariant.getSequenceName().equals(startSeqName) && leftOfMidSearchVariant.getStartPos() < startPos)) {
						//We found the start of the range since the current variant is included but the variant to left is not.
						return midSearch;
					} else {
						//The variant to the left is also part of the range we need to look further to the left
						endSearch = midSearch - 1;
					}

				}

			} else {
				//variant is to large we should search to the left
				endSearch = midSearch - 1;
			}

		}
		return variants.size();


	}

	/**
	 * Get all variant from a specific seq/chr pos to (inclusive) a specified
	 * position
	 *
	 * @param seqName
	 * @param startPos
	 * @param stopPos
	 * @return
	 */
	public Iterable<GeneticVariant> getVariantsByRange(String seqName, int startPos, int stopPos) {
		return getVariantsByRange(seqName, startPos, seqName, stopPos);
	}

	public Iterable<GeneticVariant> getVariantsBySequence(String seqName) {
		return getVariantsByRange(seqName, -1, seqName, Integer.MAX_VALUE);
	}

	/**
	 * @param seqName
	 * @param pos
	 * @return empty list if no variant found at pos
	 */
	public Iterable<GeneticVariant> getVariantAtPos(String seqName, int pos) {
		return getVariantsByRange(seqName, pos, seqName, pos);
	}

	/**
	 * Get all variants within this range. If this range is a subset of another
	 * range only the subset will be returned
	 *
	 * @return
	 */
	public List<GeneticVariant> getAllVrariantsInRange() {
		return variants;
	}

	public int size() {
		return variants.size();
	}

	@Override
	public Iterator<GeneticVariant> iterator() {
		return getAllVrariantsInRange().iterator();
	}

	public static class GeneticVariantRangeCreate {

		private final ArrayList<GeneticVariant> variants;
		private boolean isSorted;
		private GeneticVariant lastVariant = null;

		/**
		 * It is recommend to add individual variants. This constructor will
		 * always call the sort on the variants. If variants are added in order
		 * no sort will be executed
		 *
		 * @param variants
		 */
		public GeneticVariantRangeCreate(ArrayList<GeneticVariant> variants) {
			this.variants = variants;
			this.isSorted = false;
			this.lastVariant = variants.get(variants.size());
		}

		public GeneticVariantRangeCreate() {
			this.variants = new ArrayList<GeneticVariant>(20000000); // increase initial loading to match modern day datasets.
			this.isSorted = true;
		}

		public GeneticVariantRangeCreate(int expectedNumberVariants) {
			this.variants = new ArrayList<GeneticVariant>(expectedNumberVariants);
			this.isSorted = true;
		}

		public void addVariant(GeneticVariant variant) {

			if (isSorted) {
				// only check if we know that the data is actually sorted. (if it isn't. isSorted will never be changed to true)
				if (lastVariant != null && variant.getSequenceName().equals(lastVariant.getSequenceName()) && variant.getStartPos() < lastVariant.getStartPos()) {
					isSorted = false;
				} else if (lastVariant != null && ChromosomeComparator.chrASmallerChrB(variant.getSequenceName(), lastVariant.getSequenceName())) {
					isSorted = false;
				}
			}

			lastVariant = variant;
			variants.add(variant);

		}

		public GeneticVariantRange createRange() {
			if (!isSorted) {
				LOGGER.info("Variants where not sorted, loading will be faster if the data is pre-sorted");
				Collections.sort(variants, new GeneticVariantChrPosComparator());
			}
			return new GeneticVariantRange(variants);
		}

		public int size() {
			return variants.size();
		}
	}
}
