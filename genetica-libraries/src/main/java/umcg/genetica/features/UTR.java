package umcg.genetica.features;


import umcg.genetica.enums.Chromosome;

/**
 * Created by hwestra on 11/14/16.
 */
public class UTR extends Feature {

	public enum TYPE {
		PRIME3,
		PRIME5
	}

	private TYPE type;
	private Transcript transcript;

	public UTR(Chromosome chromosome, int alignmentStart, int alignmentEnd) {
		super(chromosome, alignmentStart, alignmentEnd);
	}

	public void setType(TYPE type) {
		this.type = type;
	}

	public void setTranscript(Transcript t) {
		this.transcript = t;
	}

	public TYPE getType() {
		return type;
	}

	public Transcript getTranscript() {
		return transcript;
	}
}
