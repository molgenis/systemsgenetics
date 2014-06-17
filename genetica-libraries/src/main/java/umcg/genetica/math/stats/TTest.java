/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package umcg.genetica.math.stats;

import cern.jet.random.tdouble.StudentT;
import cern.jet.stat.tdouble.Probability;

/**
 *
 * @author Harm-Jan & Marc Jan
 */
public class TTest {
    
    /**
     * Calculate t-test two-sided P-value
     * @param vals1
     * @param vals2
     * @return 
     */
    public static double test(double[] vals1, double[] vals2){
        double tTestPValue1 = -1;
        
        double mean1 = JSci.maths.ArrayMath.mean(vals1);
        double mean2 = JSci.maths.ArrayMath.mean(vals2);
        double var1 = JSci.maths.ArrayMath.variance(vals1);
        double var2 = JSci.maths.ArrayMath.variance(vals2);

        double var12 = Math.sqrt( (var1 / vals1.length) + (var2 / vals2.length) );

        double t = (mean1 - mean2) / var12;
        double df = vals1.length + vals2.length - 2;
        StudentT tDistColt = new StudentT(df, (new cern.jet.random.tdouble.engine.DRand()));

        tTestPValue1 = tDistColt.cdf(t);

        if (tTestPValue1>0.5){
            tTestPValue1 = 1 - tTestPValue1;
        }
        tTestPValue1*=2;

        return tTestPValue1;
    }
    
    /**
     * Calculate t-test Z-score
     * @param vals1
     * @param vals2
     * @return 
     */
    public static double testZscore(double[] vals1, double[] vals2){
        double tTestZScore = -1;
        
        double mean1 = JSci.maths.ArrayMath.mean(vals1);
        double mean2 = JSci.maths.ArrayMath.mean(vals2);
        double var1 = JSci.maths.ArrayMath.variance(vals1);
        double var2 = JSci.maths.ArrayMath.variance(vals2);

        double var12 = Math.sqrt( (var1 / vals1.length) + (var2 / vals2.length) );

        double t = (mean1 - mean2) / var12;
        
        double df = vals1.length + vals2.length - 2;
        StudentT tDistColt = new StudentT(df, (new cern.jet.random.tdouble.engine.DRand()));
        
        double p;
        
        if (t < 0) {
            p = tDistColt.cdf(t);
            if (p < 2.0E-323) {
                p = 2.0E-323;
            }
            tTestZScore = Probability.normalInverse(p);
        } else {
            p = tDistColt.cdf(-t);
            if (p < 2.0E-323) {
                p = 2.0E-323;
            }
            tTestZScore = -Probability.normalInverse(p);
        }

        return tTestZScore;
    }
}
