package eqtlmappingpipeline.ase;

import cern.colt.list.tint.IntArrayList;
import cern.jet.stat.tdouble.Probability;
import org.apache.log4j.Logger;

/**
 *
 * @author Patrick Deelen
 */
public class AseMle {

	private final double maxLogLikelihoodP;
	private final double maxLogLikelihood;
	private final double ratioD;
	private final double ratioP;
	private static final double[] cof = {57.1562356658629235, -59.5979603554754912,
		14.1360979747417471, -0.491913816097620199, .339946499848118887e-4,
		.465236289270485756e-4, -.983744753048795646e-4, .158088703224912494e-3,
		-.210264441724104883e-3, .217439618115212643e-3, -.164318106536763890e-3,
		.844182239838527433e-4, -.261908384015814087e-4, .368991826595316234e-5};
	private static final double[] a = new double[171];
	private static final int NTOP = 2000;
	private static final double[] aa = new double[NTOP];
	private static final Logger LOGGER = Logger.getLogger(AseMle.class);

	static {
		a[0] = 1.;
		for (int i = 1; i < 171; i++) {
			a[i] = i * a[i - 1];
		}
		for (int i = 0; i < NTOP; i++) {
			aa[i] = gammln(i + 1.);
		}
	}

	public AseMle(IntArrayList a1Counts, IntArrayList a2Counts) {

		double provisionalMaxLogLikelihood = Double.NEGATIVE_INFINITY;
		double provisionalMaxLogLikelihoodP = 0.5;
		
		
		//First calculate binominal coefficients
		double[] logBinominalCoefficients = new double[a1Counts.size()];
		for (int i = 0; i < a1Counts.size(); ++i) {
			int a1Count = a1Counts.getQuick(i);
			int totalReads = a1Count + a2Counts.getQuick(i);
			logBinominalCoefficients[i] = lnbico(totalReads, a1Count);
		}

		double logLikelihoodNull = calculateLogLikelihood(a1Counts, a2Counts, logBinominalCoefficients, 0.5);

		for (int i = 1; i <= 999; ++i) {

			double p = 0.001 * i;
			double sumLogLikelihood = calculateLogLikelihood(a1Counts, a2Counts, logBinominalCoefficients, p);

			if (sumLogLikelihood > provisionalMaxLogLikelihood) {
				provisionalMaxLogLikelihood = sumLogLikelihood;
				provisionalMaxLogLikelihoodP = p;
			}

		}

		/*
		 * Due to mathematic persion of each time adding 0.001 it might be that
		 * the null model is never really exacted tested. So now check if the
		 * null model is actually a better.
		 */
		if (logLikelihoodNull >= provisionalMaxLogLikelihood) {
			maxLogLikelihood = logLikelihoodNull;
			maxLogLikelihoodP = 0.5;
			ratioD = 0;
			ratioP = 1;
		} else {

			maxLogLikelihood = provisionalMaxLogLikelihood;
			maxLogLikelihoodP = provisionalMaxLogLikelihoodP;

			double ratioD2 = (-2d * logLikelihoodNull) + (2d * maxLogLikelihood);
			ratioD = ratioD2 < 0 ? 0 : ratioD2;
			ratioP = Probability.chiSquareComplemented(1, ratioD);

			if (Double.isInfinite(ratioD) || Double.isNaN(ratioD)) {
				LOGGER.warn("Warning invalid ratio D: " + ratioD2 + ". max log likelihood: " + maxLogLikelihood + " null log likelihood: " + logLikelihoodNull + " max log likelihood p: " + maxLogLikelihoodP);
			}

		}

	}

	private double calculateLogLikelihood(IntArrayList a1Counts, IntArrayList a2Counts, double[] logBinominalCoefficients, double p) {

		double sumLogLikelihood = 0;
		
		double logP = Math.log(p);
		double log1minP = Math.log(1 - p);

		for (int i = 0; i < a1Counts.size(); ++i) {
			sumLogLikelihood += logBinominalCoefficients[i] + (double) a1Counts.getQuick(i) * logP + (double) a2Counts.getQuick(i) * log1minP;
		}

		return sumLogLikelihood;

	}

	public double getMaxLikelihood() {
		return maxLogLikelihood;
	}

	public double getMaxLikelihoodP() {
		return maxLogLikelihoodP;
	}

	public double getRatioD() {
		return ratioD;
	}

	public double getRatioP() {
		return ratioP;
	}

	//Some borrowed functions below. Source unknown
	private static double gammln(final double xx) {
		int j;
		double x, tmp, y, ser;
		if (xx <= 0) {
			throw new IllegalArgumentException("bad arg in gammln");
		}
		y = x = xx;
		tmp = x + 5.24218750000000000; // Rational 671/128
		tmp = (x + 0.5) * Math.log(tmp) - tmp;
		ser = 0.999999999999997092;
		for (j = 0; j < 14; j++) {
			ser += cof[j] / ++y;
		}
		return tmp + Math.log(2.5066282746310005 * ser / x);
	}

	/**
	 * Returns the value n! as a floating-point number.
	 *
	 * @param n
	 * @return
	 */
	private static double factrl(final int n) {
		if (n < 0 || n > 170) {
			throw new IllegalArgumentException("factrl out of range");
		}
		return a[n];
	}

	/**
	 * Returns ln(n!)
	 *
	 * @param n
	 * @return
	 */
	private static double factln(final int n) {
		if (n < 0) {
			throw new IllegalArgumentException("negative arg in factln");
		}
		if (n < NTOP) {
			return aa[n];
		}
		return gammln(n + 1.);
	}

	/**
	 * Returns the binomial coefficient (n,k) as a floating-point number.
	 *
	 * @param n
	 * @param k
	 * @return
	 */
	protected static double bico(final int n, final int k) {
		if (n < 0 || k < 0 || k > n) {
			throw new IllegalArgumentException("bad args in bico");
		}
		if (n < 171) {
			return Math.floor(0.5 + factrl(n) / (factrl(k) * factrl(n - k)));
		}
		return Math.floor(0.5 + Math.exp(factln(n) - factln(k) - factln(n - k)));
	}

	protected static double lnbico(final int n, final int k) {
		if (n < 0 || k < 0 || k > n) {
			throw new IllegalArgumentException("bad args in bico");
		}
		if (n < 171) {
			return Math.log(Math.floor(0.5 + factrl(n) / (factrl(k) * factrl(n - k))));
		}
		return factln(n) - factln(k) - factln(n - k);
	}
}
