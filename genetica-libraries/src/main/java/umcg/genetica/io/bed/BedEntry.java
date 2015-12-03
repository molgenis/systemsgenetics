package umcg.genetica.io.bed;

import umcg.genetica.variantAnnotator.GenomicRange;

/**
 *
 * @author Patrick Deelen
 */
public class BedEntry implements GenomicRange {

	private final String chr;
	private final int startPosition;
	private final int stopPositions;
	private final String name;
	private final double score;

	public BedEntry(String chr, int startPosition, int stopPositions) {
		this.chr = chr;
		this.startPosition = startPosition;
		this.stopPositions = stopPositions;
		this.name = null;
		this.score = Double.NaN;
	}

	public BedEntry(String chr, int startPosition, int stopPositions, String name, double score) {
		this.chr = chr;
		this.startPosition = startPosition;
		this.stopPositions = stopPositions;
		this.name = name;
		this.score = score;
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

	@Override
	public String getSeqname() {
		return chr;
	}

	public String getName() {
		return name;
	}

	public double getScore() {
		return score;
	}
	
	

}
