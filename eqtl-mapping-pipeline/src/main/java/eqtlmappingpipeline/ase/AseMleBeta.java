package eqtlmappingpipeline.ase;

import cern.colt.list.tint.IntArrayList;
import cern.jet.stat.tdouble.Probability;
import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import org.apache.log4j.Logger;

/**
 *
 * @author Patrick Deelen
 */
public class AseMleBeta {

	private final double maxLogLikelihoodP;
	private final double maxLogLikelihood;
	private final double maxLogLikelihoodTheta;
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
	private static final Logger LOGGER = Logger.getLogger(AseMleBeta.class);
	private static final BigDecimal probabilityStep = new BigDecimal(0.001);
	private static final MathContext mathContextProbability = new MathContext(3, RoundingMode.HALF_UP);//Must fit step percision
	protected static final double[] probabilities;
	protected static final double[] oneMinProbabilities;
	private static final int maxThetaToTest = 1000;

	static {
		a[0] = 1.;
		for (int i = 1; i < 171; i++) {
			a[i] = i * a[i - 1];
		}
		for (int i = 0; i < NTOP; i++) {
			aa[i] = gammln(i + 1.);
		}
	}

	static {
		BigDecimal one = new BigDecimal(1);

		int stepCount = ((int) (1 / probabilityStep.doubleValue())) - 1;

		probabilities = new double[stepCount];
		oneMinProbabilities = new double[stepCount];

		int i = 0;
		for (BigDecimal p = probabilityStep; p.compareTo(one) < 0; p = p.add(probabilityStep)) {

			double p2 = p.round(mathContextProbability).doubleValue();
			probabilities[i] = p2;
			oneMinProbabilities[i] = 1 - p2;

			++i;
		}

	}

	public AseMleBeta(IntArrayList a1Counts, IntArrayList a2Counts) {

		double provisionalMaxLogLikelihood = Double.NEGATIVE_INFINITY;
		double provisionalMaxLogLikelihoodP = 0.5;
		double provisionalMaxLogLikelihoodTheta = Double.NaN;


		//First calculate binominal coefficients
		double[] logBinominalCoefficients = new double[a1Counts.size()];
		for (int i = 0; i < a1Counts.size(); ++i) {
			int a1Count = a1Counts.getQuick(i);
			int totalReads = a1Count + a2Counts.getQuick(i);
			logBinominalCoefficients[i] = lnbico(totalReads, a1Count);
		}

		for (int i = 0; i < probabilities.length; ++i) {

			for (int theta = 1; theta <= maxThetaToTest; ++theta) {

				double sumLogLikelihood = 0;
				for (int s = 0; s < a1Counts.size(); ++s) {
					double alfa = theta * oneMinProbabilities[i];
					double beta = theta * probabilities[i];
					sumLogLikelihood += logBinominalCoefficients[s] + betaln(a2Counts.getQuick(s) + alfa, a1Counts.getQuick(s) + beta) - betaln(alfa, beta);
				}

				if (sumLogLikelihood > provisionalMaxLogLikelihood) {
					provisionalMaxLogLikelihood = sumLogLikelihood;
					provisionalMaxLogLikelihoodP = probabilities[i];
					provisionalMaxLogLikelihoodTheta = theta;
				}

			}

		}

		double logLikelihoodNullTheta = provisionalMaxLogLikelihoodTheta;
		double logLikelihoodNull = 0;

		for (int s = 0; s < a1Counts.size(); ++s) {
			double alfa = logLikelihoodNullTheta * 0.5;
			double beta = logLikelihoodNullTheta * 0.5;
			logLikelihoodNull += logBinominalCoefficients[s] + betaln(a2Counts.getQuick(s) + alfa, a1Counts.getQuick(s) + beta) - betaln(alfa, beta);
		}

		if (Double.isInfinite(logLikelihoodNull)) {
			throw new RuntimeException("Something went wrong during ASE analysis. This should nog happen, please contact developers");
		}

		//Make sure to use null model in case of tie
		if (logLikelihoodNull >= provisionalMaxLogLikelihood) {

			if (provisionalMaxLogLikelihoodTheta == maxThetaToTest) {
				LOGGER.warn("Warning, possibly did not reach optimum for theta in beta binominal model");
			}

			maxLogLikelihood = logLikelihoodNull;
			maxLogLikelihoodP = 0.5;
			maxLogLikelihoodTheta = logLikelihoodNullTheta;
			ratioD = 0;
			ratioP = 1;
		} else {

			maxLogLikelihood = provisionalMaxLogLikelihood;
			maxLogLikelihoodP = provisionalMaxLogLikelihoodP;
			maxLogLikelihoodTheta = provisionalMaxLogLikelihoodTheta;

			if (maxLogLikelihoodTheta == maxThetaToTest) {
				LOGGER.warn("Warning, possibly did not reach optimum for theta in beta binominal model");
			}

			double ratioD2 = (-2d * logLikelihoodNull) + (2d * maxLogLikelihood);
			ratioD = ratioD2 < 0 ? 0 : ratioD2;
			ratioP = Probability.chiSquareComplemented(1, ratioD);

			if (Double.isInfinite(ratioD) || Double.isNaN(ratioD)) {
				LOGGER.warn("Warning invalid ratio D: " + ratioD2 + ". max log likelihood: " + maxLogLikelihood + " null log likelihood: " + logLikelihoodNull + " max log likelihood p: " + maxLogLikelihoodP);
			}

		}

	}

	public double getMaxLikelihood() {
		return maxLogLikelihood;
	}

	public double getMaxLikelihoodP() {
		return maxLogLikelihoodP;
	}

	public double getMaxLogLikelihoodTheta() {
		return maxLogLikelihoodTheta;
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

	//based on numerical recipies
	private static double betaln(final double z, final double w) {
		return gammln(z) + gammln(w) - gammln(z + w);
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
