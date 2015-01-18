package umcg.genetica.io.binInteraction;

/**
 *
 * @author Patrick Deelen
 */
public class BinaryInteractionQtlZscores {
	
	private final double[] zscore;
	private final int[] sampleCounts;
	private final double metaZscore;
	
	private final static double[] emptyDoubleArray = new double[0];
	private final static int[] emptyIntArray = new int[0];

	public BinaryInteractionQtlZscores(double[] zscore, int[] sampleCounts, double metaZscore) {
		this.zscore = zscore;
		this.sampleCounts = sampleCounts;
		this.metaZscore = metaZscore;
	}
	
	public BinaryInteractionQtlZscores(double[] zscore, int[] sampleCounts) {
		this.zscore = zscore;
		this.sampleCounts = sampleCounts;
		this.metaZscore = Double.NaN;
	}
	
	public BinaryInteractionQtlZscores(double metaZscore) {
		this.zscore = emptyDoubleArray;
		this.sampleCounts = emptyIntArray;
		this.metaZscore = metaZscore;
	}

	public double[] getZscore() {
		return zscore;
	}

	public int[] getSampleCounts() {
		return sampleCounts;
	}

	public double getMetaZscore() {
		return metaZscore;
	}

}
