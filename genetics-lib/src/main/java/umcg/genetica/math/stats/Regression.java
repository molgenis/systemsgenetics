/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats;

/**
 *
 * @author harmjan
 */
public class Regression {

    public static double[] getLinearRegressionCoefficients(double[] x, double[] y) {
        if (y.length != x.length) {
            throw new IllegalStateException("Error in linear regression: X and Y should have the same length! " + x.length + " vs " + y.length);
        }
        double meanX = Descriptives.mean(x);
        double meanY = Descriptives.mean(y);
        double sXY = 0;
        double sXX = 0;
        double sYY = 0;
        for (int i = 0; i < x.length; i++) {
            double xmin = x[i] - meanX;
            double ymin = y[i] - meanY;
            double xminsqrd = xmin * xmin;
            double yminsqrd = ymin * ymin;
            sXY += (xmin * ymin);
            sXX += xminsqrd;
            sYY += yminsqrd;
        }

        double beta = sXY / sXX;
        double alpha = meanY - beta * meanX;

        double ssxy = 0;
        for (int i = 0; i < y.length; i++) {
            double yexp = alpha + (beta * x[i]);
            ssxy += ((y[i] - yexp) * (y[i] - yexp));
        }

//        double se = sYY - ((sXY * sXY) / (sXX * sXX));
//        
//        
//        
//        se /= (x.length - 2);
//        se = Math.sqrt(se);

        double se = Math.sqrt((ssxy / (x.length - 2))) / Math.sqrt(sXX);
        double t = beta / se;
        return new double[]{beta, alpha, se, t};
    }
}
