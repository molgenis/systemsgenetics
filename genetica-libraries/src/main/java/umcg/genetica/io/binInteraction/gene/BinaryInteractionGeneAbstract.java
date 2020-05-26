package umcg.genetica.io.binInteraction.gene;

/**
 *
 * @author Patrick Deelen
 */
public abstract class BinaryInteractionGeneAbstract implements BinaryInteractionGene {

	private final String name;
	private final String chr;
	private final int start;
	private final int end;

	public BinaryInteractionGeneAbstract(String name, String chr, int start, int end) {
		this.name = name;
		this.chr = chr;
		this.start = start;
		this.end = end;
	}

	@Override
	public String getName() {
		return name;
	}

	@Override
	public String getChr() {
		return chr;
	}

	@Override
	public int getStart() {
		return start;
	}

	@Override
	public int getEnd() {
		return end;
	}
	
	
	
}
