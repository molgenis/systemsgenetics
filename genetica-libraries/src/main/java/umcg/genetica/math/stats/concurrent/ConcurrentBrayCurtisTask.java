/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats.concurrent;

import java.util.concurrent.Callable;
import umcg.genetica.containers.Pair;
import umcg.genetica.math.stats.Descriptives;

/**
 *
 * @author harmjan
 */
public class ConcurrentBrayCurtisTask implements Callable<Pair<Integer, double[]>> {

    private final int x;
    private final double[][] matrix1;

    // correlate rows on matrix1
    public ConcurrentBrayCurtisTask(double[][] matrix1, int x) {
        this.matrix1 = matrix1;
        this.x = x;
    }

    // WARNING: this method only calculates half of the matrix.
    @Override
    public Pair<Integer, double[]> call() throws Exception {
        // determine what to correlate
        
        double[] results = new double[matrix1.length];
        double[] xarr = matrix1[x];
        for (int i = x+1; i < matrix1.length; i++) {
            double[] yarr = matrix1[i];
            double r = Descriptives.BrayCurtisDistance(xarr, yarr);
            results[i] = r;
        }
        results[x] = 0;
        return new Pair<Integer, double[]>(x, results);
    }
}
