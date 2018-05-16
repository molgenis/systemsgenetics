/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.exception.NoDataException;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.stat.ranking.NaNStrategy;
import org.apache.commons.math3.stat.ranking.NaturalRanking;
import org.apache.commons.math3.stat.ranking.TiesStrategy;
import org.apache.commons.math3.util.FastMath;

/**
 * Adapted version of org.apache.commons.math3.stat.inference.MannWhitneyUTest
 *
 *
 * @author patri
 */
public class MannWhitneyUTest2 {

	double u1 = Double.NaN;
	double u2 = Double.NaN;
	double u = Double.NaN;
	double p = Double.NaN;
	boolean isPset = false;
	double z = Double.NaN;
	int n1 = 0;
	int n2 = 0;
	double auc = 0;

	/**
	 * Ranking algorithm.
	 */
	private final NaturalRanking naturalRanking;

	public MannWhitneyUTest2() {
		naturalRanking = new NaturalRanking(NaNStrategy.FIXED,
				TiesStrategy.AVERAGE);
	}

	/**
	 * Create a test instance using the given strategies for NaN's and ties.
	 * Only use this if you are sure of what you are doing.
	 *
	 * @param nanStrategy specifies the strategy that should be used for
	 * Double.NaN's
	 * @param tiesStrategy specifies the strategy that should be used for ties
	 */
	public MannWhitneyUTest2(final NaNStrategy nanStrategy,
			final TiesStrategy tiesStrategy) {
		naturalRanking = new NaturalRanking(nanStrategy, tiesStrategy);
	}

	/**
	 * Ensures that the provided arrays fulfills the assumptions.
	 *
	 * @param x first sample
	 * @param y second sample
	 * @throws NullArgumentException if {@code x} or {@code y} are {@code null}.
	 * @throws NoDataException if {@code x} or {@code y} are zero-length.
	 */
	private void ensureDataConformance(final double[] x, final double[] y)
			throws NullArgumentException, NoDataException {

		if (x == null
				|| y == null) {
			throw new NullArgumentException();
		}
		if (x.length == 0
				|| y.length == 0) {
			throw new NoDataException();
		}
	}

	/**
	 * Concatenate the samples into one array.
	 *
	 * @param x first sample
	 * @param y second sample
	 * @return concatenated array
	 */
	private double[] concatenateSamples(final double[] x, final double[] y) {
		final double[] z = new double[x.length + y.length];

		System.arraycopy(x, 0, z, 0, x.length);
		System.arraycopy(y, 0, z, x.length, y.length);

		return z;
	}

	/**
	 * Computes the <a
	 * href="http://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U"> Mann-Whitney
	 * U statistic</a> comparing mean for two independent samples possibly of
	 * different length.
	 * <p>
	 * This statistic can be used to perform a Mann-Whitney U test evaluating
	 * the null hypothesis that the two independent samples has equal mean.
	 * </p>
	 * <p>
	 * Let X<sub>i</sub> denote the i'th individual of the first sample and
	 * Y<sub>j</sub> the j'th individual in the second sample. Note that the
	 * samples would often have different length.
	 * </p>
	 * <p>
	 * <strong>Preconditions</strong>:
	 * <ul>
	 * <li>All observations in the two samples are independent.</li>
	 * <li>The observations are at least ordinal (continuous are also
	 * ordinal).</li>
	 * </ul>
	 * </p>
	 *
	 * @param x the first sample
	 * @param y the second sample
	 * @throws NullArgumentException if {@code x} or {@code y} are {@code null}.
	 * @throws NoDataException if {@code x} or {@code y} are zero-length.
	 */
	public void setData(final double[] x, final double[] y)
			throws NullArgumentException, NoDataException {

		isPset = false;
		p = Double.NaN;

		ensureDataConformance(x, y);

		final double[] xy = concatenateSamples(x, y);
		final double[] ranks = naturalRanking.rank(xy);

		double sumRankX = 0;

		n1 = x.length;
		n2 = y.length;

		/*
         * The ranks for x is in the first x.length entries in ranks because x
         * is in the first x.length entries in z
		 */
		for (int i = 0; i < n1; ++i) {
			sumRankX += ranks[i];
		}

		/*
         * U1 = R1 - (n1 * (n1 + 1)) / 2 where R1 is sum of ranks for sample 1,
         * e.g. x, n1 is the number of observations in sample 1.
		 */
		u1 = sumRankX - ((long) n1 * (n1 + 1)) / 2;

		final long n1n2prod = (long) n1 * n2;

		/*
         * It can be shown that U1 + U2 = n1 * n2
		 */
		u2 = n1n2prod - u1;

		u = FastMath.max(u1, u2);

		// http://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U#Normal_approximation
		final double EU = n1n2prod / 2.0;
		final double VarU = n1n2prod * (n1 + n2 + 1) / 12.0;

		z = (u2 - EU) / FastMath.sqrt(VarU);

		auc = u1 / n1n2prod;

	}

	public double getU1() {
		return u1;
	}

	public double getU2() {
		return u2;
	}

	public double getP() {
		if (!isPset) {
			// No try-catch or advertised exception because args are valid
//			final NormalDistribution standardNormal = new NormalDistribution(0, 1);

//			p = 2 * standardNormal.cumulativeProbability(-FastMath.abs(z));
			double x = FastMath.abs(z);
			double[] b = {0.319381530, -0.356563782, 1.781477937, -1.821255978, 1.330274429};
			double p = 0.2316419;
			double t = 1 / (1 + p * x);
			double fact = t;
			double sum = 0;
			for (int i = 0; i <= b.length - 1; i++) {
				sum += b[i] * fact;
				fact *= t;
			}
			p = 2d * sum * Math.exp(-x * x / 2.0d) / (Math.sqrt(2d * Math.PI));
			this.p = p;

		}
		return p;
	}

	public double getZ() {
		return z;
	}

	public NaturalRanking getNaturalRanking() {
		return naturalRanking;
	}

	public int getN1() {
		return n1;
	}

	public int getN2() {
		return n2;
	}

	public double getU() {
		return u;
	}

	public double getAuc() {
		return auc;
	}

}
