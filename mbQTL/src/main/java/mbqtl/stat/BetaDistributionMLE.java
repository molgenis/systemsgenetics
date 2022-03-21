package mbqtl.stat;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import org.apache.commons.math3.special.Beta;

public class BetaDistributionMLE {


	private final static double BETA_SHAPE1_MIN = 0.1;
	private final static double BETA_SHAPE2_MIN = 1;
	private final static double BETA_SHAPE1_MAX = 10;
	private final static double BETA_SHAPE2_MAX = 1000000;


	public class BetaDistributionMLEFunction implements MultivariateFunction {
		private double[] p;

		public void setP(double[] p) {
			this.p = p;
		}

		@Override
		public double value(double[] doubles) {
			double beta_shape1 = doubles[0];
			double beta_shape2 = doubles[1];
			double diff = -1.0 * ((beta_shape1 - 1) * p[0] + (beta_shape2 - 1) * p[1] - p[2] * Beta.logBeta(beta_shape1, beta_shape2));
			// System.out.println(beta_shape1 + "\t" + beta_shape2 + "\t" + diff);
			return diff;
		}
	}


	public double[] fit(double[] pvals) {
		double mean = JSci.maths.ArrayMath.mean(pvals);
		double var = JSci.maths.ArrayMath.variance(pvals);
		double beta_shape1 = mean * (mean * (1 - mean) / var - 1);
		double beta_shape2 = beta_shape1 * (1 / mean - 1);
		return fit(pvals, beta_shape1, beta_shape2);
	}

	public double[] fit(double[] pvals, double alpha, double beta) {

//		System.out.println();
//		System.out.println("Fit beta: ");
//		System.out.println("Start alpha: " + alpha);
//		System.out.println("Start beta: " + beta);
		double[] x = new double[]{alpha, beta};
		double[] ss = new double[]{alpha / 10, beta / 10};
		double[] par = new double[3];
		for (int q = 0; q < pvals.length; q++) {
			if (pvals[q] == 1) {
				pvals[q] = 0.999999;
			}
			par[0] += Math.log(pvals[q]);
			par[1] += Math.log(1 - pvals[q]);
		}
		par[2] = pvals.length;

		BetaDistributionMLEFunction func = new BetaDistributionMLEFunction();
		func.setP(par);

		double relthresh = 0.0001;
		double absthresh = 0.0001;
		SimplexOptimizer optimizer = new SimplexOptimizer(relthresh, absthresh);
		final PointValuePair optimum = optimizer.optimize(
				new MaxEval(1000),
				new ObjectiveFunction(func),
				GoalType.MINIMIZE,
				new InitialGuess(x),
				new NelderMeadSimplex(ss));

		return optimum.getPoint();
	}


}
