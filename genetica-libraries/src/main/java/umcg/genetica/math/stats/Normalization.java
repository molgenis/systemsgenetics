package umcg.genetica.math.stats;

/**
 *
 * @author juha
 */
public class Normalization {

    private Normalization() {
    }

    public static double[] standardNormalize(double[] vals) {
        
        double[] normalized = new double[vals.length];
        double mean = Descriptives.mean(vals);
        double stdev = Descriptives.variance(vals, mean);
        for (int i = 0; i < vals.length; i++) {
            normalized[i] -= mean;
            normalized[i] /= stdev;
        }
        return normalized;
    }
}
