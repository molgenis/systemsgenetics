package umcg.genetica.io.binInteraction;

/**
 *
 * @author Patrick Deelen
 */
public class BinaryInteractionCohort {
	
	private final String name;
	private final int sampleCount;

	public BinaryInteractionCohort(String name, int sampleCount) {
		this.name = name;
		this.sampleCount = sampleCount;
	}

	public String getName() {
		return name;
	}

	public int getSampleCount() {
		return sampleCount;
	}

}
