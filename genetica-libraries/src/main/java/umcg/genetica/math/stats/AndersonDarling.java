package umcg.genetica.math.stats;

import JSci.maths.statistics.NormalDistribution;
import java.util.logging.Level;
import java.util.logging.Logger;
import umontreal.iro.lecuyer.probdist.AndersonDarlingDist;

/**
 *
 * Anderson-Darling test for normality. Uses SSJ (GPL) for A-D distribution.
 * 
 * http://en.wikipedia.org/wiki/Anderson%E2%80%93Darling_test
 * ref 10th Aug 2012
 * 
 * @author juha
 * 
 * @warning. The ssj dependency, depended on colt_1.2.0. This dependency is no longer there. If in the background colt functions are used this could lead to errors.
 * 
 * 
 */
public class AndersonDarling {

    private static final Logger LOGGER = Logger.getLogger(AndersonDarling.class.getName());
    private static final NormalDistribution norm = new NormalDistribution();

    public static enum IS_DATA_NORMALIZED {

        YES, NO
    };

    // no instantiation, all static methods
    private AndersonDarling() {
    }

    public static double getASquared(double[] vals) {
        return getASquared(vals, IS_DATA_NORMALIZED.NO);
    }

    public static double getASquared(double[] vals, IS_DATA_NORMALIZED idn) {

        if (vals == null || vals.length < 3) {
            throw new IllegalArgumentException("An array of at least 3 values should be given.");
        }

        if (idn == IS_DATA_NORMALIZED.NO) {
            LOGGER.log(Level.INFO, "Standard normalizing data.");
            vals = Normalization.standardNormalize(vals);
            LOGGER.log(Level.INFO, "Data normalized.");
        }

        double n = vals.length;
        double a2 = 0; // A^2 test statistic
        for (int i = 1; i <= vals.length; i++) {
            a2 += (2 * i - 1) * Math.log(norm.cumulative(vals[i - 1])) + (2 * (n - i) + 1) * Math.log(1 - norm.cumulative(vals[i - 1]));
        }
        a2 = -n - a2 / n;

        return a2;
    }

    public static double getP(double[] vals, IS_DATA_NORMALIZED idn) {
        double asq = getASquared(vals, idn);
        return getP(asq, vals.length);
    }

    public static double getP(double asq, int n) {
        return AndersonDarlingDist.cdf(n, asq);
    }
}
