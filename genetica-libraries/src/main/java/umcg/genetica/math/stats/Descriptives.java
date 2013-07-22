/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats;

/**
 *
 * @author harmjan
 */
public class Descriptives {

    public static double[] m_zScoreToPValue;
    public static double[] m_sqrtSample;

    /**
     * Calculate the mean for an array of floats
     *
     * @param v The array of floats to calculate the mean for
     * @return The mean of the float array
     */
    public static double mean(float[] v) {
        double sum = 0;
        for (int k = 0; k < v.length; k++) {
            sum += v[k];
        }
        return (sum / (double) v.length);
    }

    /**
     * Calculate the variance for an array of floats
     *
     * @param v The array of floats for which the variance should be calculated
     * @param mean The mean of the array of floats
     * @return The variance for the array of floats
     */
    public static double variance(float[] v, double mean) {
        double ans = 0.0;
        for (int i = 0; i < v.length; i++) {
            ans += (v[i] - mean) * (v[i] - mean);
        }
        return ans / (v.length - 1);
    }

    public static void lookupSqrt(int totalNrSamples) {
        //Fast look-up service for squared root of number of samples:
        double[] sqrtSample = new double[totalNrSamples + 1];
        for (int nrSamples = 0; nrSamples <= totalNrSamples; nrSamples++) {
            sqrtSample[nrSamples] = Math.sqrt((double) nrSamples);
        }
        m_sqrtSample = sqrtSample;
    }

    public static void zScoreToPValue() {
        //Fast look-up service for normal distribution:
        JSci.maths.statistics.NormalDistribution normDist = new JSci.maths.statistics.NormalDistribution();
        double[] zScoreToPValue = new double[300001];
        for (int z = 0; z <= 300000; z++) {
            double zScore = ((double) z - 150000d) / 5000d;
            double pValue = 0;
            //Avoid complementation:
            if (zScore > 0) {
                pValue = normDist.cumulative(-zScore);
            } else {
                pValue = normDist.cumulative(zScore);
            }
            //Two sided test:
            if (pValue > 0.5) {
                pValue = 1 - pValue;
            }
            pValue *= 2.0d;
            //Eventual P Value is stored:
            zScoreToPValue[z] = pValue;
        }
        m_zScoreToPValue = zScoreToPValue;
    }

    public static double zScore(double value, double mean, double variance) {
        if (variance > 0.0 && mean > 0.0) {
            double sd = Math.sqrt(variance);
            double val = value - mean;
            double zscore = val / sd;
            zscore = Math.sqrt(Math.pow(zscore, 2));
            return zscore;
        } else {
            return Double.MAX_VALUE;
        }
    }

    public static double mean(double[] v) {
        double sum = 0;
        for (int k = 0; k < v.length; k++) {
            sum += v[k];
        }
        return (sum / (double) v.length);
    }

    public static double variance(double[] v) {
        double mean = mean(v);
        return variance(v, mean);
    }

    public static double variance(double[] v, double mean) {
        double ans = 0.0;
        for (int i = 0; i < v.length; i++) {
            ans += (v[i] - mean) * (v[i] - mean);
        }
        return ans / (v.length - 1);
    }

    public static int getZScorePvalueIndex(double zScore) {
        int zScoreIndex = (int) (zScore * 5000.0d) + 150000;
        if (zScoreIndex < 0) {
            zScoreIndex = 0;
        }
        if (zScoreIndex > 300000) {
            zScoreIndex = 300000;
        }
        return zScoreIndex;
    }

    public static double convertZscoreToPvalue(double zScore) {

        int zScoreIndex = (int) (zScore * 5000.0d) + 150000;
        if (zScoreIndex < 0) {
            zScoreIndex = 0;
        }
        if (zScoreIndex > 300000) {
            zScoreIndex = 300000;
        }

        if (m_zScoreToPValue == null) {
            zScoreToPValue();
        }
        double pValueOverall = m_zScoreToPValue[zScoreIndex];
        return pValueOverall;
    }

    public static double getSqrt(int nrTotalSamples) {
        if (m_sqrtSample == null) {
            System.out.println("ERROR: square-root table not correctly initialized.");
            System.exit(-1);
        }
        return m_sqrtSample[nrTotalSamples];
    }

    public static double sum(double[] v) {
        double sum = 0;
        for (double d : v) {
            sum += d;
        }
        return sum;
    }

    public static double absSum(double[] v) {
        double sum = 0;
        for (double d : v) {
            sum += Math.abs(d);
        }
        return sum;
    }
}
