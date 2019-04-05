/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats;

import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.util.stream.IntStream;

/**
 * @author harmjan
 */
public class Log2Transform {

    public static void log2transform(DoubleMatrixDataset<String, String> rawData) {
        double minValue = Double.MAX_VALUE;

        int probeCount = rawData.rows();
        int sampleCount = rawData.columns();

        for (int p = 0; p < probeCount; p++) {
            for (int s = 0; s < sampleCount; s++) {
                double v = rawData.getElementQuick(p, s);
                if (v < minValue) {
                    minValue = v;
                }
            }
        }

        System.out.println("\nLog2 transforming data: Absolute minimum Expression Value:\t" + minValue);
        double multiplier = 1.0d / Math.log10(2.0d);
        double finalMinValue = minValue;
        IntStream.range(0, probeCount).parallel().forEach(p -> {
            for (int s = 0; s < sampleCount; s++) {
                if (finalMinValue <= 0) {
                    double v = Math.log10(rawData.getElementQuick(p, s) - finalMinValue + 1) * multiplier;
                    rawData.setElementQuick(p, s, v);
                } else {
                    double v = Math.log10(rawData.getElementQuick(p, s) + 1) * multiplier;
                    rawData.setElementQuick(p, s, v);
                }
            }
        });
    }
}
