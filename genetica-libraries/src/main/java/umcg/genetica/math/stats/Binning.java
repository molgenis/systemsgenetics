/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats;

/**
 *
 * @author juha
 */
public class Binning {

    public static int[] getBins(double[] data, double min, double max, int nrBins) {

        int[] bins = new int[nrBins];
        double step = (max - min) / nrBins;
        for (double d : data) {
            if (d < min) {
                continue;
            }
            for (int b = 0; b < nrBins; b++) {
                if (d < min + (b + 1) * step) {
                    bins[b]++;
                    break;
                }
            }
        }
        return bins;
    }
}
