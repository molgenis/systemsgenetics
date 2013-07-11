/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats;

/**
 *
 * @author harmjan
 */
public class KolmogorovSmirnov {

    public static double getProbabilityKolmogorovSmirnov(double d, double n1, double n2) { //Numerical recipes code

        double ne = (n1 * n2) / (n1 + n1);
        double alam = ((Math.sqrt(ne) + 0.12 + 0.11 / Math.sqrt(ne)) * d);
        double eps1 = 0.001;
        double eps2 = 1E-8;
        double a2 = -2 * alam * alam;
        double fac = 2;
        double probKS = 0;
        double termbf = 0;
        for (int j = 1; j <= 100; j++) {
            double term = fac * Math.exp(a2 * (double) j * (double) j);
            probKS += term;
            if (Math.abs(term) <= eps1 * termbf || Math.abs(term) <= eps2 * probKS) {
                return probKS;
            }
            fac = -fac;
            termbf = Math.abs(term);
        }
        return 1;
    }

    /**
     * 
     * Calculates a 2-sample Kolmogorov-Smirnov p-value.
     * 
     * @param vals1 sample 1, must be sorted!
     * @param vals2 sample 2, must be sorted!
     * @return K-S p
     */
    public static double getKSP(double[] vals1, double[] vals2) {
        double z = getKSZ(vals1, vals2);
        return zToP(z);
    }

    /**
     * 
     * Calculates a 2-sample Kolmogorov-Smirnov Z value. Note! Arrays must be sorted.
     * 
     * @param vals1 sample 1, must be sorted.
     * @param vals2 sample 2, must be sorted.
     * @return K-S Z
     */
    public static double getKSZ(double[] vals1, double[] vals2) {

        int j1 = 0, j2 = 0;
        double d = 0.0, d1, d2, dt, en1, en2, fn1 = 0.0, fn2 = 0.0;

        en1 = vals1.length;
        en2 = vals2.length;

        while (j1 < en1 && j2 < en2) {
            if (vals1[j1] == Double.NaN || vals2[j2] == Double.NaN) {
                throw new IllegalArgumentException("No NaNs allowed!");
            }
            if ((d1 = vals1[j1]) <= (d2 = vals2[j2])) {
                fn1 = ++j1 / en1;
            }

            if (d2 <= d1) {
                fn2 = ++j2 / en2;
            }

            if ((dt = Math.abs(fn2 - fn1)) > d) {
                d = dt;
            }
        }

        double z = Math.sqrt(en1 * en2 / (en1 + en2)) * d;
        return z;
    }

    /**
     * 
     * Converts a Kolmogorov-Smirnov Z value to a p-value.
     * 
     * @param z
     * @return p
     */
    public static double zToP(double z) {

        double EPS1 = 1.0e-6, EPS2 = 1.0e-16;
        int j;
        double a2, fac = 2.0, sum = 0.0, term, termbf = 0.0;

        a2 = -2.0 * z * z;

        for (j = 1; j <= 100; ++j) {
            term = fac * Math.exp(a2 * j * j);
            sum += term;
            if (Math.abs(term) <= EPS1 * termbf || Math.abs(term) <= EPS2 * sum) {
                return sum;
            }
            fac = -fac;
            termbf = Math.abs(term);
        }

        return 1.0;
    }
}
