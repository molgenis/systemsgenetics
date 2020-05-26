/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats.concurrent;

import java.util.concurrent.Callable;
import umcg.genetica.containers.Pair;
import umcg.genetica.math.stats.Correlation;

/**
 *
 * @author marc jan
 */
public class ConcurrentCovariationTask implements Callable<Pair<Integer, double[]>> {

    private final int x;
    private final double[] meanOfSamples;
    private final double[][] matrix1;

    // correlate rows on matrix1
    public ConcurrentCovariationTask(double[][] matrix1, double[] meanOfSamples, int x) {
        this.matrix1 = matrix1;
        this.meanOfSamples = meanOfSamples;
        this.x = x;
    }

    // WARNING: this method only calculates half of the matrix.
    @Override
    public Pair<Integer, double[]> call() throws Exception {
        // determine what to covariate
        
        double[] results = new double[matrix1.length];
        double[] xarr = matrix1[x];
        for (int i = x; i < matrix1.length; i++) {
            double[] yarr = matrix1[i];
            double r = Correlation.covariate(meanOfSamples[x], meanOfSamples[i], xarr, yarr);
            results[i] = r;
        }
//        results[x] = 1;
        return new Pair<Integer, double[]>(x, results);
    }
}
