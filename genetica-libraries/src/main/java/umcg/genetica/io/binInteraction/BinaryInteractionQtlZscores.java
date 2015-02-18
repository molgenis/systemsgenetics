package umcg.genetica.io.binInteraction;

/**
 *
 * @author Patrick Deelen
 */
public class BinaryInteractionQtlZscores {
	
	private final double[] zscores;
	private final int[] sampleCounts;
	private final double metaZscore;
	
	private final static double[] emptyDoubleArray = new double[0];
	private final static int[] emptyIntArray = new int[0];

	public BinaryInteractionQtlZscores(double[] zscores, int[] sampleCounts, double metaZscore) throws BinaryInteractionFileException {
		this.zscores = zscores;
		this.sampleCounts = sampleCounts;
		this.metaZscore = metaZscore;
		checkArrays();
	}
	
	public BinaryInteractionQtlZscores(double[] zscores, int[] sampleCounts) throws BinaryInteractionFileException {
		this.zscores = zscores;
		this.sampleCounts = sampleCounts;
		this.metaZscore = Double.NaN;
		checkArrays();
	}
	
	public BinaryInteractionQtlZscores(double zscore, int sampleCounts) throws BinaryInteractionFileException {
		this.zscores = new double[] {zscore};
		this.sampleCounts = new int[] {sampleCounts};
		this.metaZscore = Double.NaN;
		checkArrays();
	}
	
	public BinaryInteractionQtlZscores(double metaZscore) {
		this.zscores = emptyDoubleArray;
		this.sampleCounts = emptyIntArray;
		this.metaZscore = metaZscore;
	}
	
	private void checkArrays() throws BinaryInteractionFileException{
				
		if(zscores.length != sampleCounts.length){
			throw new BinaryInteractionFileException("All arrays for qtl z-scores must be equal length");
		}
		
	}

	public double[] getZscores() {
		return zscores;
	}

	public int[] getSampleCounts() {
		return sampleCounts;
	}

	public double getMetaZscore() {
		return metaZscore;
	}
	
	public int getCohortCount(){
		return zscores.length;
	}

}
