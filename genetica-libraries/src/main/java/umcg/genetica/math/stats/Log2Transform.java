/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats;

import java.util.stream.IntStream;

/**
 *
 * @author harmjan
 */
public class Log2Transform {
    
    public static void log2transform(double[][] rawData){
        double minValue = Double.MAX_VALUE;
        
        int probeCount = rawData.length;
        int sampleCount = rawData[probeCount-1].length;
        
        for (int p=0; p<probeCount; p++) {
            for (int s=0; s<sampleCount; s++) {
                if (rawData[p][s]<minValue) {
                    minValue = rawData[p][s];
                }
            }
        }
        
        System.out.println("\nLog2 transforming data: Absolute minimum Expression Value:\t" + minValue);
        double multiplier = 1.0d / Math.log10(2.0d);
        double finalMinValue = minValue;
        IntStream.range(0,probeCount).parallel().forEach(p->{
            for (int s=0; s<sampleCount; s++) {
                if (finalMinValue <= 0) {
                    rawData[p][s] = (double) (Math.log10(rawData[p][s] - finalMinValue + 1) * multiplier);
                } else {
                    rawData[p][s] = (double) (Math.log10(rawData[p][s] + 1) * multiplier);
                }
            }
        });
    }
}
