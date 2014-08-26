package umcg.genetica.io.bedgraph;

import umcg.genetica.variantAnnotator.GenomicRange;

/**
 *
 * @author Patrick Deelen
 */
public class BedGraphEntry implements GenomicRange {

	private final String chr;
	private final int startPosition;
	private final int stopPositions;
	private final double value;

	public BedGraphEntry(String chr, int startPosition, int stopPositions, double value) {
		this.chr = chr;
		this.startPosition = startPosition;
		this.stopPositions = stopPositions;
		this.value = value;
	}

	@Override
	public int getStart() {
		return this.startPosition;
	}

	@Override
	public int getEnd() {
		return this.stopPositions - 1;
	}

	/**
	 *
	 * @return original bed file end of range.
	 */
	public int getEndExclusive() {
		return this.stopPositions;
	}

	public double getValue() {
		return value;
	}

	@Override
	public String getSeqname() {
		return chr;
	}

}
