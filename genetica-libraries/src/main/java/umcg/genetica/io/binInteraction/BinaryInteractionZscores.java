package umcg.genetica.io.binInteraction;

import java.util.Arrays;

/**
 *
 * @author Patrick Deelen
 */
public class BinaryInteractionZscores {

	private final int[] samplesInteractionCohort;
	private final double[] zscoreSnpCohort;
	private final double[] zscoreCovariateCohort;
	private final double[] zscoreInteractionCohort;
	private final double[] rSquaredCohort;
	private final double[] zscoreInteractionFlippedCohort;
	private final double zscoreSnpMeta;
	private final double zscoreCovariateMeta;
	private final double zscoreInteractionMeta;
	private final double zscoreInteractionFlippedMeta;
	
	protected final static double[] emptyDoubleArray = new double[0];
	protected final static int[] emptyIntArray = new int[0];

	public BinaryInteractionZscores(int[] samplesInteractionCohort, double[] zscoreSnpCohort, double[] zscoreCovariateCohort, double[] zscoreInteractionCohort, double[] rSquaredCohort, double[] zscoreInteractionFlippedCohort, double zscoreSnpMeta, double zscoreCovariateMeta, double zscoreInteractionMeta, double zscoreInteractionFlippedMeta) throws BinaryInteractionFileException {
		this.samplesInteractionCohort = samplesInteractionCohort;
		this.zscoreSnpCohort = zscoreSnpCohort;
		this.zscoreCovariateCohort = zscoreCovariateCohort;
		this.zscoreInteractionCohort = zscoreInteractionCohort;
		this.rSquaredCohort = rSquaredCohort;
		this.zscoreInteractionFlippedCohort = zscoreInteractionFlippedCohort;
		this.zscoreSnpMeta = zscoreSnpMeta;
		this.zscoreCovariateMeta = zscoreCovariateMeta;
		this.zscoreInteractionMeta = zscoreInteractionMeta;
		this.zscoreInteractionFlippedMeta = zscoreInteractionFlippedMeta;
		checkArrays();
	}	
	
	public BinaryInteractionZscores(int[] samplesInteractionCohort, double[] zscoreSnpCohort, double[] zscoreCovariateCohort, double[] zscoreInteractionCohort, double[] rSquaredCohort, double[] zscoreInteractionFlippedCohort) throws BinaryInteractionFileException {
		this.samplesInteractionCohort = samplesInteractionCohort;
		this.zscoreSnpCohort = zscoreSnpCohort;
		this.zscoreCovariateCohort = zscoreCovariateCohort;
		this.zscoreInteractionCohort = zscoreInteractionCohort;
		this.rSquaredCohort = rSquaredCohort;
		this.zscoreInteractionFlippedCohort = zscoreInteractionFlippedCohort;
		this.zscoreSnpMeta = Double.NaN;
		this.zscoreCovariateMeta = Double.NaN;
		this.zscoreInteractionMeta = Double.NaN;
		this.zscoreInteractionFlippedMeta = Double.NaN;
		checkArrays();
	}
	
	public BinaryInteractionZscores(int[] samplesInteractionCohort, double[] zscoreSnpCohort, double[] zscoreCovariateCohort, double[] zscoreInteractionCohort, double[] rSquaredCohort) throws BinaryInteractionFileException {
		this.samplesInteractionCohort = samplesInteractionCohort;
		this.zscoreSnpCohort = zscoreSnpCohort;
		this.zscoreCovariateCohort = zscoreCovariateCohort;
		this.zscoreInteractionCohort = zscoreInteractionCohort;
		this.rSquaredCohort = rSquaredCohort;
		this.zscoreInteractionFlippedCohort = new double[samplesInteractionCohort.length];
		Arrays.fill(zscoreInteractionFlippedCohort, Double.NaN);
		this.zscoreSnpMeta = Double.NaN;
		this.zscoreCovariateMeta = Double.NaN;
		this.zscoreInteractionMeta = Double.NaN;
		this.zscoreInteractionFlippedMeta = Double.NaN;
		checkArrays();
	}
	
	public BinaryInteractionZscores(int samplesInteractionCohort, double zscoreSnpCohort, double zscoreCovariateCohort, double zscoreInteractionCohort, double rSquaredCohort, double zscoreInteractionFlippedCohort) throws BinaryInteractionFileException {
		this.samplesInteractionCohort = new int[] {samplesInteractionCohort};
		this.zscoreSnpCohort = new double[] {zscoreSnpCohort};
		this.zscoreCovariateCohort = new double[] {zscoreCovariateCohort};
		this.zscoreInteractionCohort = new double[] {zscoreInteractionCohort};
		this.rSquaredCohort = new double[] {rSquaredCohort};
		this.zscoreInteractionFlippedCohort = new double[] {zscoreInteractionFlippedCohort};
		this.zscoreSnpMeta = Double.NaN;
		this.zscoreCovariateMeta = Double.NaN;
		this.zscoreInteractionMeta = Double.NaN;
		this.zscoreInteractionFlippedMeta = Double.NaN;
		checkArrays();
	}

	public BinaryInteractionZscores(int[] samplesInteractionCohort, double[] zscoreSnpCohort, double[] zscoreCovariateCohort, double[] zscoreInteractionCohort, double[] rSquaredCohort, double zscoreSnpMeta, double zscoreCovariateMeta, double zscoreInteractionMeta) throws BinaryInteractionFileException {
		this.samplesInteractionCohort = samplesInteractionCohort;
		this.zscoreSnpCohort = zscoreSnpCohort;
		this.zscoreCovariateCohort = zscoreCovariateCohort;
		this.zscoreInteractionCohort = zscoreInteractionCohort;
		this.rSquaredCohort = rSquaredCohort;
		this.zscoreInteractionFlippedCohort = new double[samplesInteractionCohort.length];
		Arrays.fill(zscoreInteractionFlippedCohort, Double.NaN);
		this.zscoreSnpMeta = zscoreSnpMeta;
		this.zscoreCovariateMeta = zscoreCovariateMeta;
		this.zscoreInteractionMeta = zscoreInteractionMeta;
		this.zscoreInteractionFlippedMeta = Double.NaN;
		checkArrays();
	}
	
	public BinaryInteractionZscores(double zscoreSnpMeta, double zscoreCovariateMeta, double zscoreInteractionMeta, double zscoreInteractionFlippedMeta) {
		this.samplesInteractionCohort = emptyIntArray;
		this.zscoreSnpCohort = emptyDoubleArray;
		this.zscoreCovariateCohort = emptyDoubleArray;
		this.zscoreInteractionCohort = emptyDoubleArray;
		this.rSquaredCohort = emptyDoubleArray;
		this.zscoreInteractionFlippedCohort = emptyDoubleArray;
		this.zscoreSnpMeta = zscoreSnpMeta;
		this.zscoreCovariateMeta = zscoreCovariateMeta;
		this.zscoreInteractionMeta = zscoreInteractionMeta;
		this.zscoreInteractionFlippedMeta = zscoreInteractionFlippedMeta;
	}
	
	public BinaryInteractionZscores(double zscoreSnpMeta, double zscoreCovariateMeta, double zscoreInteractionMeta) {
		this.samplesInteractionCohort = emptyIntArray;
		this.zscoreSnpCohort = emptyDoubleArray;
		this.zscoreCovariateCohort = emptyDoubleArray;
		this.zscoreInteractionCohort = emptyDoubleArray;
		this.rSquaredCohort = emptyDoubleArray;
		this.zscoreInteractionFlippedCohort = emptyDoubleArray;
		this.zscoreSnpMeta = zscoreSnpMeta;
		this.zscoreCovariateMeta = zscoreCovariateMeta;
		this.zscoreInteractionMeta = zscoreInteractionMeta;
		this.zscoreInteractionFlippedMeta = Double.NaN;
	}
	
	private void checkArrays() throws BinaryInteractionFileException{
		
		int count = samplesInteractionCohort.length;
		
		if(
				zscoreSnpCohort.length != count
				|| zscoreCovariateCohort.length != count
				|| zscoreInteractionCohort.length != count
				|| rSquaredCohort.length != count
				|| zscoreInteractionFlippedCohort.length != count){
			throw new BinaryInteractionFileException("All arrays for interaction z-scores must be equal length");
		}
		
	}

	public int[] getSamplesInteractionCohort() {
		return samplesInteractionCohort;
	}

	public double[] getZscoreSnpCohort() {
		return zscoreSnpCohort;
	}

	public double[] getZscoreCovariateCohort() {
		return zscoreCovariateCohort;
	}

	public double[] getZscoreInteractionCohort() {
		return zscoreInteractionCohort;
	}

	public double[] getrSquaredCohort() {
		return rSquaredCohort;
	}

	public double[] getZscoreInteractionFlippedCohort() {
		return zscoreInteractionFlippedCohort;
	}

	public double getZscoreSnpMeta() {
		return zscoreSnpMeta;
	}

	public double getZscoreCovariateMeta() {
		return zscoreCovariateMeta;
	}

	public double getZscoreInteractionMeta() {
		return zscoreInteractionMeta;
	}

	public double getZscoreInteractionFlippedMeta() {
		return zscoreInteractionFlippedMeta;
	}
	
	public int getCohortCount(){
		return zscoreInteractionCohort.length;
	}
	
}
