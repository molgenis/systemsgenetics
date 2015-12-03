package nl.systemsgenetics.eqtlinteractionanalyser.eqtlinteractionanalyser;

/**
 *
 * @author Patrick Deelen
 */
public class GeneAnnotation {

	private final String ensg;
	private final String huho;
	private final String chr;
	private final int start;
	private final int end;

	public GeneAnnotation(String ensg, String huho, String chr, int start, int end) {
		this.ensg = ensg;
		this.huho = huho;
		this.chr = chr;
		this.start = start;
		this.end = end;
	}

	public String getEnsg() {
		return ensg;
	}

	public String getHuho() {
		return huho;
	}

	public String getChr() {
		return chr;
	}

	public int getStart() {
		return start;
	}

	public int getEnd() {
		return end;
	}
	
}
