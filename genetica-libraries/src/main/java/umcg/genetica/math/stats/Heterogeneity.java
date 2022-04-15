package umcg.genetica.math.stats;

import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.util.Primitives;

import java.util.ArrayList;

/**
 * @author harmjan
 */
public class Heterogeneity {

    public static Triple<Double, Double, Integer> getISq(ArrayList<Double> isqd, ArrayList<Integer> isqn) {

        double[] z = Primitives.toPrimitiveArr(isqd.toArray(new Double[0]));
        int[] n = Primitives.toPrimitiveArr(isqn.toArray(new Integer[0]));

        return getISq(z, n);
    }

    /*
     * This class was freely adapted from Abecasis' METAL package for meta-analysis
     * (http://www.sph.umich.edu/csg/abecasis/Metal/)
     * it should produce output that is highly correlated with the isqFromCorrAndNFixedEffect function below
     * datasetZ: array of z-scores per dataset, NaN for missing
     * datasetWeights: array of sample sizes (not squared)
     */
    public static Triple<Double, Double, Integer> getISq(double[] datasetZ, int[] datasetWeights) {

        double weightedZ = 0;
        int totalSample = 0;
        for (int d = 0; d < datasetZ.length; d++) {
            if (!Double.isNaN(datasetZ[d])) {
                weightedZ += Math.sqrt(datasetWeights[d]) * datasetZ[d];
                totalSample += datasetWeights[d];
            }
        }

        double hetSum = 0;
        int hetDf = 0;
        for (int d = 0; d < datasetZ.length; d++) {
            if (!Double.isNaN(datasetZ[d])) {
                double expectedZ = Math.sqrt(datasetWeights[d]) * weightedZ / totalSample;

                hetSum += (datasetZ[d] - expectedZ) * (datasetZ[d] - expectedZ);
                hetDf++;
            }
        }

        double p = 1d;
        double i = 0d;

        double iExp = ((hetSum - hetDf + 1) / hetSum) * 100d;
        if (hetDf <= 1 || hetSum < 1E-7) {
            p = 1;
        } else {
            p = ChiSquare.getP(hetDf - 1, hetSum);
        }


        if (hetSum <= (hetDf - 1) || hetDf <= 1) {
            i = 0;
        } else {
            i = iExp;
        }

        if (i > 0) {
            i /= 100;
        }
        return new Triple<Double, Double, Integer>(i, p, totalSample);
    }


    // Borenstein et al, Introduction to Meta-analysis, 2009
    // https://onlinelibrary.wiley.com/doi/book/10.1002/9780470743386
    // pp 98++
    // r: array of correlation coefficients
    // n: array of sample sizes
    public static double[] isqFromCorrAndNFixedEffect(double[] r, int[] n) {
        double sumW = 0;
        double sumWY = 0;
        double sumWY2 = 0;
        double sumW2 = 0;
        int nrDs = 0;
        for (int d = 0; d < r.length; d++) {
            if (!Double.isNaN(r[d])) {
                double y = 0.5 * Math.log((1 + r[d]) / (1 - r[d])); // fisher z-transform of r to effect size
                double var = 1d / (n[d] - 3);                       // estimated variance of effect size
                double se = Math.sqrt(var);
                double w = 1d / var;
                double wy = w * y;
                double wy2 = (w * (y * y));
                double w2 = w * w;

                sumW += w;
                sumWY += wy;
                sumWY2 += wy2; // incorrect
                sumW2 += w2;
                nrDs++;
            }
        }

        double M = sumWY / sumW;                      // meta-analysis effect size
        double varM = 1d / sumW;                      // meta-analysis effect size variance
        double seM = Math.sqrt(varM);
        double lowerCI = M - 1.96 * seM;
        double upperCI = M + 1.96 * seM;
        double Z = M / seM;                           // meta-analysis Z-score
        double Q = sumWY2 - ((sumWY * sumWY) / sumW); // cochran's Q
        double df = nrDs - 1;
        double C = sumW - (sumW2 / sumW);
        double T2 = (Q - df) / C;                     // tau-squared
        double I2 = ((Q - df) / Q) * 100;             // heterogeneity

        // I2 can be negative for casaes of extremely low heterogeneity
        if (I2 < 0) {
            I2 = 0;
        }

        // see page 124 for confidence intervals around I2
        double rM = (Math.exp(2 * M) - 1) / (Math.exp(2 * M) + 1); // inverse Fisher- z transform; meta-analysed correlation coefficient
        double rMLL = (Math.exp(2 * lowerCI) - 1) / (Math.exp(2 * lowerCI) + 1);  // Lower CI for rM
        double rMUL = (Math.exp(2 * upperCI) - 1) / (Math.exp(2 * upperCI) + 1);  // Upper CI for rM
        double or = Math.exp(M); // OR conversion

        return new double[]{M, varM, seM, Z, Q, C, T2, I2};
    }
}
