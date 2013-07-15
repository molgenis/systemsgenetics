/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats;

import JSci.maths.ArrayMath;
import cern.colt.Arrays;
import java.io.IOException;
import org.apache.commons.collections.primitives.ArrayDoubleList;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author juha
 */
public class TheilSen {

    private TheilSen() {
    }
    
    public static void main(String[] args) throws IOException {
        TextFile tf = new TextFile("/Users/juha/Downloads/correlationTest.txt", false);
        int lines = tf.countLines() - 1;
        String line = tf.readLine();
        double[] a1 = new double[lines];
        double[] a2 = new double[lines];
        int i = 0;
        while ((line = tf.readLine()) != null) {
            String[] split = line.split("\t");
            a1[i] = Double.parseDouble(split[1]);
            a2[i] = Double.parseDouble(split[2]);
            i++;
        }
        double[] descriptives = getDescriptives(a1, a2);
        System.out.println(Arrays.toString(descriptives));
        double[] linearRegressionCoefficients = Regression.getLinearRegressionCoefficients(a1, a2);
        System.out.println(Arrays.toString(linearRegressionCoefficients));
    }

    /**
     * 
     * Calculates Theil-Sen estimation for given arrays.
     * 
     * @param v1
     * @param v2
     * @return [y intercept, slope, 2.5 percentile, 97.5 percentile, number of pairs]
     */
    public static double[] getDescriptives(double[] v1, double[] v2) {

        if (v1.length != v2.length) {
            throw new IllegalArgumentException("Arrays must be of the same length! " + v1.length + ", " + v2.length);
        }

        ArrayDoubleList slopesList = new ArrayDoubleList();
        int cnt = 0;
        for (int i = 0; i < v1.length; i++) {
            double x = v1[i];
            double y = v2[i];
            for (int j = i + 1; j < v1.length; j++) {
                if (x != v1[j]) { // x must be different, otherwise slope becomes infinite
                    double slope = (v2[j] - y) / (v1[j] - x);
                    slopesList.add(slope);
                    ++cnt;
                }
            }
        }

        double[] slopes = slopesList.toArray();
        double median1 = ArrayMath.median(v1);
        double median2 = ArrayMath.median(v2);
        double slope = ArrayMath.median(slopes);
        double yI = median2 - slope * median1;
        double p1 = ArrayMath.percentile(slopes, 0.025d);
        double p2 = ArrayMath.percentile(slopes, 0.975d);

        return new double[]{yI, slope, p1, p2, cnt};

    }
}
