/*
 * This class was freely translated from R BioConductor package "qvalue"
 * Storey JD. (2002) A direct approach to false discovery rates. Journal of the Royal Statistical Society, Series B, 64: 479-498.
 * Storey JD and Tibshirani R. (2003) Statistical significance for genome-wide studies. Proceedings of the National Academy of Sciences, 100: 9440-9445.
 * Storey JD. (2003) The positive false discovery rate: A Bayesian interpretation and the q-value. Annals of Statistics, 31: 2013-2035.
 * Storey JD, Taylor JE, and Siegmund D. (2004) Strong control, conservative point estimation, and simultaneous conservative consistency of false discovery rates: A unified approach. Journal of the Royal Statistical Society, Series B, 66: 187-205.
 */
package umcg.genetica.math.stats;

import java.util.ArrayList;
import java.util.Arrays;
import umcg.genetica.math.interpolation.CubicSpline;
import umcg.genetica.util.Primitives;
import umcg.genetica.util.RankDoubleArray;

/**
 *
 * @author harmjan
 */
public class QValues {

    public static enum SMOOTHER {

        BOOTSTRAP, SMOOTHER
    };

    public static double[] qvalue(double[] pvals, Double fdr) {
        double[] lambda = initlambda();
        return qvalue(pvals, lambda, fdr, false, 3, SMOOTHER.SMOOTHER, false);
    }

    public static double[] qvalue(double[] pvals, double[] lambda, Double fdr, boolean robust, Integer smootfhDf, SMOOTHER smoother, boolean smoothLog) {

        if (pvals == null || pvals.length == 0) {
            throw new IllegalArgumentException("Qvalue - Error: list of p-values is null or length equals zero");
        }

        double maxp = Primitives.max(pvals);
        double minp = Primitives.min(pvals);
        if (maxp > 1 || minp < 0) {
            throw new IllegalArgumentException("Qvalue - Error: pvals should be between 0 and 1. Max: " + maxp + "\tMin: " + minp);
        }

        if (lambda == null) {
            lambda = initlambda();
        } else if (lambda.length > 1 && lambda.length < 4) {
            throw new IllegalArgumentException("Qvalue - Error: if length of lambda > 1, lambda should have at least 4 values");
        } else if (lambda.length > 1 && (Primitives.max(lambda) > 1 || Primitives.min(lambda) < 0)) {
            throw new IllegalArgumentException("Qvalue - Error: lambda should be between [0,1)");
        } else if (lambda.length == 1 && (lambda[0] < 0 || lambda[0] > 1)) {
            throw new IllegalArgumentException("Qvalue - Error: lambda should be between [0,1)");
        }

        int m = pvals.length;
        double pi0 = 0;

        // pi0 estimation
        if (lambda.length == 1) {


            double lval = lambda[0];
            double nrvalsltlambda = 0;
            for (double p : pvals) {
                if (p >= lval) {
                    nrvalsltlambda += 1d;
                }
            }

            nrvalsltlambda /= pvals.length;
            pi0 = nrvalsltlambda / (1 - lval);
            pi0 = Primitives.min(new double[]{pi0, 1d});

        } else {
            double[] pi0Intermediate = new double[lambda.length];
            for (int i = 0; i < pi0Intermediate.length; i++) {
                double lval = lambda[i];
                double nrvalsltlambda = 0;
                for (double p : pvals) {
                    if (p >= lval) {
                        nrvalsltlambda += 1d;
                    }
                }
                pi0Intermediate[i] = (nrvalsltlambda / lambda.length) / (1 - lval);
            }
            if (smoother == SMOOTHER.SMOOTHER) {

                if (smoothLog) {
                    for (int i = 0; i < pi0Intermediate.length; i++) {
                        pi0Intermediate[i] = Math.log(pi0Intermediate[i]); // natural log
                    }
                }

                CubicSpline spline = new CubicSpline(lambda, pi0Intermediate);
                pi0 = spline.interpolate(Primitives.max(lambda));

                if (smoothLog) {
                    pi0 = Math.exp(pi0);
                }

                pi0 = Primitives.min(new double[]{pi0, 1d});

            } else if (smoother == SMOOTHER.BOOTSTRAP) {
                double minpi0 = Primitives.min(pi0Intermediate);
                double[] mse = new double[lambda.length];

                double[] piboot = new double[lambda.length];
                for (int i = 0; i < 100; i++) {

                    // sample m values from pvals, put them in piboot. (or actually randomize??)
                    double[] pboot = sampleWithReplacement(pvals, pvals.length);
                    for (int j = 0; j < lambda.length; j++) {
                        double lval = lambda[i];

                        double nrvalsltlambda = 0;
                        for (double p : pboot) {
                            if (p > lval) {
                                nrvalsltlambda += 1d;
                            }
                        }

                        nrvalsltlambda /= pboot.length;
                        piboot[i] = nrvalsltlambda / (1 - lval);
                    }
                    for (int j = 0; j < mse.length; j++) {
                        double d = (piboot[j] - minpi0);
                        mse[j] += (d * d);
                    }
                }

                double minmse = Primitives.min(mse);

                ArrayList<Integer> minMSEIndices = new ArrayList<Integer>();
                for (int j = 0; j < mse.length; j++) {
                    if (mse[j] == minmse) {
                        minMSEIndices.add(j);
                    }
                }
                double[] minpi0s = new double[minMSEIndices.size()];
                for (int j = 0; j < minpi0s.length; j++) {
                    minpi0s[j] = pi0Intermediate[minMSEIndices.get(j)];
                }
                pi0 = Primitives.min(minpi0s);
                pi0 = Primitives.min(new double[]{pi0, 1d});

            } else {
                throw new IllegalArgumentException("SMOOTHER argument must be either be BOOTSTRAP or SMOOTHER");
            }



        }

        // now check whether pi0 was calculated correctly
        if (pi0 <= 0) {
            System.err.println("ERROR: The estimated pi0 <= 0. Check that you have valid p-values or use another lambda estimator");
        }

        if (fdr != null && (fdr <= 0 || fdr > 1)) {
            throw new IllegalArgumentException("ERROR: FDR level should be between 0 and 1");
        }

        RankDoubleArray rda = new RankDoubleArray();
        double[] rankedpvals = rda.rank(pvals);

        double[] v = specialQValueRank(pvals);

        

        return null;
    }

    private static double[] initlambda() {
        double[] lambda = new double[18];
        for (int i = 0; i < lambda.length; i++) {
            lambda[i] = (double) i * 0.05;
        }
        return lambda;
    }

    private static double[] sampleWithReplacement(double[] vals, int numToSample) {
        double[] output = new double[numToSample];
        for (int i = 0; i < numToSample; i++) {
            int randomPos = (int) Math.floor(Math.random() * vals.length);
            output[i] = vals[randomPos];
        }
        return output;
    }

    private static double[] specialQValueRank(double[] pvals) {
        throw new UnsupportedOperationException("Not yet implemented");
    }
}
