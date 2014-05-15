/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.methylation;

import cern.colt.matrix.tdouble.DoubleMatrix2D;

/**
 *
 * @author MarcJan
 */
public class ConvertBetaAndMvalues {
    
    public static void transformToMvalue(double[][] rawData){
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
    
    public static void transformMToBetavalue(DoubleMatrix2D rawData){
        
        int probeCount = rawData.rows();
        int sampleCount = rawData.columns();
        
        for (int p=0; p<probeCount; p++) {
            for (int s=0; s<sampleCount; s++) {
                double tmpBeta = Math.pow(2, rawData.getQuick(p, s));
                rawData.setQuick(p, s, tmpBeta/(tmpBeta+1));
            }
        }
    }
    
    public static void transformMToBetavalue(double[][] rawData){
        
        int probeCount = rawData.length;
        int sampleCount = rawData[probeCount-1].length;
        
        for (int p=0; p<probeCount; p++) {
            for (int s=0; s<sampleCount; s++) {
                double tmpBeta = Math.pow(2, rawData[p][s]);
                rawData[p][s] = tmpBeta/(tmpBeta+1);
            }
        }
    }
    
}
