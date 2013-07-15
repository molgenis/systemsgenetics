/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.methylation;

/**
 *
 * @author MarcJan
 */
public class ConvertBetaToMvalue {
    
    public static void transToMvalue(double[][] rawData){
        double minValue = Double.MAX_VALUE;
        double maxValue = Double.MIN_VALUE;
        
        int probeCount = rawData.length;
        int sampleCount = rawData[probeCount-1].length;
        
        for (int p=0; p<probeCount; p++) {
            for (int s=0; s<sampleCount; s++) {
                if (rawData[p][s]!=0  && rawData[p][s]<minValue) {
                    minValue = rawData[p][s];
                }
                if (rawData[p][s]!=1  && rawData[p][s]>maxValue) {
                    maxValue = rawData[p][s];
                }
            }
        }
        
        double multiplier = 1.0d / Math.log10(2.0d);
        minValue = (double) (Math.log10((minValue/(1-minValue))) * multiplier);
        maxValue = (double) (Math.log10((maxValue/(1-maxValue))) * multiplier);
        
        for (int p=0; p<probeCount; p++) {
            for (int s=0; s<sampleCount; s++) {
                if(rawData[p][s] == 0){
                    rawData[p][s] = minValue;
                } else if(rawData[p][s] == 1){
                    rawData[p][s] = maxValue;
                } else {
                    rawData[p][s] = (double) (Math.log10((rawData[p][s] / (1 - rawData[p][s]))) * multiplier);
                }
            }
        }
    }
}
