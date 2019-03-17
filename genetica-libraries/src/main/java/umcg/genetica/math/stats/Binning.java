/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats;

/**
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

    public static double[] getFreq(int[] cts, int total) {
        double[] output = new double[cts.length];
        for (int d = 0; d < output.length; d++) {
            if (cts[d] == 0) {
                output[d] = 0;
            } else {
                output[d] = (double) cts[d] / total;
            }
        }
        return output;
    }

    public static double[] getBounds(double min, double max, int nrbins) {
        double stepsize = (max - min) / nrbins;
        double[] output = new double[nrbins];
        for (int d = 0; d < output.length; d++) {
            if (d == 0) {
                output[d] = min;
            } else {
                output[d] = output[d - 1] + stepsize;
            }
        }
        return output;
    }

    public static double[] getBoundsBinCenter(double min, double max, int nrbins) {
        double stepsize = (max - min) / nrbins;
        double[] output = new double[nrbins];
        for (int d = 0; d < output.length; d++) {
            if (d == 0) {
                output[d] = min + (stepsize/2);
            } else {
                output[d] = output[d - 1] + stepsize;
            }
        }
        return output;
    }
}
