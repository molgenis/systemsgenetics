package mbqtl.stat;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;

public class DFMLE {

    public class DegreeOfFreedomFunction implements MultivariateFunction {
        private double[] corr;

        public void setData(double[] corr) {
            this.corr = corr;
        }

        @Override
        public double value(double[] doubles) {

            double[] pvals = new double[corr.length];
            double mean = 0;
            for (int c = 0; c < corr.length; c++) {
                pvals[c] = PVal.getPvalue(corr[c], Math.abs(doubles[0])); // don't allow guess to be negative.
                mean += pvals[c];
            }
            mean /= pvals.length;
            double variance = 0.0;
            for (int c = 0; c < pvals.length; c++) {
                variance += (pvals[c] - mean) * (pvals[c] - mean);
            }
            variance /= (pvals.length - 1);

            double shape2 = Math.abs((mean * (mean * (1 - mean) / variance - 1)) - 1.0);
//            System.out.println(doubles[0] + "\t" + mean + "\t" + variance + "\t" + shape2);
            //cout << "O = " << mean << " " << shape2 << endl;
            return shape2;
        }
    }

    public Double learnDF(double[] bestPermCorrelations, double trueDF) {
        DegreeOfFreedomFunction func = new DegreeOfFreedomFunction();
        func.setData(bestPermCorrelations);
        double[] x = new double[]{trueDF};
        double[] ss = new double[]{trueDF * 0.1};

        // provide bounds to keep simplex solution positive
//        MultivariateFunctionMappingAdapter adapter = new MultivariateFunctionMappingAdapter(func, new double[]{ss[0]}, new double[]{trueDF * 3});

        double relthresh = 0.01;
        double absthresh = 0.01;
        SimplexOptimizer optimizer = new SimplexOptimizer(relthresh, absthresh);
        final PointValuePair optimum = optimizer.optimize(
                new MaxEval(20),
                new ObjectiveFunction(func),
                GoalType.MINIMIZE,
                new InitialGuess(x),
                new NelderMeadSimplex(ss));

//        final PointValuePair optimum = optimizer.optimize(
//                new MaxEval(20),
//                new ObjectiveFunction(func),
//                GoalType.MINIMIZE,
//                new InitialGuess(ss),
//                new NelderMeadSimplex(x));


        return Math.abs(optimum.getPoint()[0]);
    }


}
