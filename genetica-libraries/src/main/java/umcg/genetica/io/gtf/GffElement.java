/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.gtf;

import java.util.Collections;
import java.util.Map;
import umcg.genetica.variantAnnotator.GenomicRange;

/**
 *
 * @author Patrick Deelen
 */
public class GffElement implements GenomicRange{
	
	//http://mblab.wustl.edu/GTF22.html
	
	/**
	 * name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix.
	 */
	private final String seqname;
	
	/**
	 * name of the program that generated this feature, or the data source (database or project name)
	 */
	private final String source;
	
	/**
	 * feature type name, e.g. Gene, Variation, Similarity
	 */
	private final String feature;
	
	/**
	 * Start position of the feature, with sequence numbering starting at 1
	 */
	private final int start;
	
	/**
	 * End position of the feature, with sequence numbering starting at 1
	 */
	private final int end;
	
	/**
	 * A floating point value
	 */
	private final float score;
	
	/**
	 * defined as + (forward) or - (reverse)
	 */
	private final char strand;
	
	/**
	 * One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on
	 */
	private final int frame;
	
	/**
	 * A semicolon-separated list of tag-value pairs, providing additional information about each feature
	 */
	private final Map<String, String> attributes;

	public GffElement(String seqname, String source, String feature, int start, int end, float score, char strand, int frame, Map<String, String> attributes) {
		this.seqname = seqname;
		this.source = source;
		this.feature = feature;
		this.start = start;
		this.end = end;
		this.score = score;
		this.strand = strand;
		this.frame = frame;
		this.attributes = Collections.unmodifiableMap(attributes);
	}

	@Override
	public String getSeqname() {
		return seqname;
	}

	public String getSource() {
		return source;
	}

	public String getFeature() {
		return feature;
	}

	@Override
	public int getStart() {
		return start;
	}

	@Override
	public int getEnd() {
		return end;
	}

	public float getScore() {
		return score;
	}

	public char getStrand() {
		return strand;
	}

	public int getFrame() {
		return frame;
	}

	public Map<String, String> getAttributes() {
		return attributes;
	}
	
	public boolean hasAttribute(String attributeName){
		return attributes.containsKey(attributeName);
	}
	
	public String getAttributeValue(String attributeName){
		return attributes.get(attributeName);
	}

	@Override
	public String toString() {
		StringBuilder builder = new StringBuilder();
		builder.append("Seq: ");
		builder.append(getSeqname());
		builder.append("\nPos: ");
		builder.append(getStart());
		builder.append(" to ");
		builder.append(getEnd());
		builder.append("\nsource: ");
		builder.append(getSource());
		builder.append("\nfeature: ");
		builder.append(getFeature());
		builder.append("\nscore: ");
		builder.append(getScore());
		builder.append("\nstrand: ");
		builder.append(getStrand());
		builder.append("\nframe: ");
		builder.append(getFrame());
		builder.append("\nattributes: ");
		builder.append(attributes.toString());
		return builder.toString();
		
	}
	
}
