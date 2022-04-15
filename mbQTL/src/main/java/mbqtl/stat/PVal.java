package mbqtl.stat;

import JSci.maths.statistics.NormalDistribution;
import cern.jet.random.tdouble.StudentT;
import cern.jet.random.tdouble.engine.DRand;
import org.apache.commons.math3.distribution.FDistribution;

public class PVal {
    private static final NormalDistribution NORM_DIST = new NormalDistribution(0.0D, 1.0D);

    public static double getPvalueFDist(double corr, double df) {
        FDistribution f = new FDistribution(1, df);
        double fstat = df * corr * corr / (1 - corr * corr);
        return 1 - f.cumulativeProbability(fstat);
    }

    public static double getPvalue(double spearman, double df) {
        if (Math.abs(spearman - -1.0D) < 1.0E-7D) {
            spearman = -0.9999D;
        }

        if (Math.abs(spearman - 1.0D) < 1.0E-7D) {
            spearman = 0.9999D;
        }
        double t = spearman / Math.sqrt((1.0D - spearman * spearman) / (double) (df));
        StudentT tDistColt = new StudentT((double) (df), new DRand());
        double p = 1;
        if (spearman >= 0) {
            p = tDistColt.cdf(-t) * 2;
        } else {
            p = tDistColt.cdf(t) * 2;
        }
        if (p < 2.0E-323D) {
            p = 2.0E-323D;
        }
        return p;
    }

    public static double zToP(double z) {
        double p;
        if (z > 0.0D) {
            p = NORM_DIST.cumulative(-z);
        } else {
            p = NORM_DIST.cumulative(z);
        }

        if (p > 0.5D) {
            p = 1.0D - p;
        }

        p *= 2.0D;
        if (p < 2.0E-323D) {
            p = 2.0E-323D;
        }
        return p;
    }


}
