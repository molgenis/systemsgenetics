/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package eqtlmappingpipeline.mixupmapper.stat.dist;

import eqtlmappingpipeline.graphics.histogram.FrequencyDistributionHistogram;
import eqtlmappingpipeline.graphics.histogram.Histogram;
import java.util.ArrayList;
import java.util.Collections;



/**
 *
 * @author harm-jan
 */
public class DiscreteDist extends Dist {
    private double  binInterval;
    private int     numTotalValues;
    private Bin[]   bins;
    private int     binIterator;
    
    public DiscreteDist(){
        super();
        binInterval     = 0.0;
        numTotalValues  = 0;
        binIterator     = 0;
    }

   public void createDist(double _min, double _max, int numBins){
        min     = _min;
        max     = _max;
        createBins(numBins);
    }

    public void createDist(ArrayList<Double> alValues, int numBins){
        double _max     = Collections.max(alValues);
        double _min     = Collections.min(alValues);

        createDist(alValues, numBins, _min, _max);
    }

    public void createDist(ArrayList<Double> alValues, int numBins, double _min, double _max){
        min     = _min;
        max     = _max;

        
        createBins(numBins);

        if(alValues != null){
            // put values in bins
            double  dBinNum = 0.0;
            int     iBinNum = 0;
            for (Object value: alValues){
                double dValue = (Double) value;
                dBinNum = (dValue - min) / binInterval;
                
                iBinNum = (int) Math.floor(dBinNum);

                if(iBinNum >= bins.length) {
                    iBinNum--;
                }
                if(iBinNum < bins.length){
                    bins[iBinNum].add();
                    numTotalValues++;
                } else {
                    System.out.println("could not include value: "+value+" because it falls out of the range: "+min+"-"+max);
                }
            }
            // calcStatistics(alValues);
        }
    }

    public double getBinInterval(){
        return binInterval;
    }

    public int getTotalNumValues(){
        return numTotalValues;
    }

    public void setTotalNumValues(int total){
        numTotalValues = total;
    }

    public int getNumBins(){
        return bins.length;
    }

    public void calcCumulative(){
        int add    = 0;
        for(int binNum =0; binNum< bins.length; binNum++){
            if(bins[binNum] != null){
                bins[binNum].setCumulative(0);
                add += bins[binNum].getCount();
                bins[binNum].setCumulative(add);
            }
        }
    }

    public void plot(String outfile){
        Histogram hist = new Histogram(512,512);
        hist.plot(this);
        hist.draw(outfile);
    }

    public void plotFreqDist(String outfile){
        FrequencyDistributionHistogram hist = new FrequencyDistributionHistogram(512,512);
        hist.plot(this);
        hist.draw(outfile);
    }

    public boolean hasNext(){
        if(binIterator < bins.length){
            return true;
        } else {
            return false;
        }
    }

    public Bin getNext(){
        if(binIterator < bins.length){
            Bin currentBin = bins[binIterator];
            binIterator++;
            return currentBin;
        } else {
            return null;
        }
    }

     public Bin getPrev(){
        if(binIterator - 1 >= 0 && binIterator - 1 < bins.length-1){
            Bin currentBin = bins[binIterator - 1];
            return currentBin;
        } else {
            return null;
        }
    }

    public void resetIterator(){
        binIterator = 0;
    }

    public void convertToFrequency(int total){
        for(int binNum =0; binNum< bins.length; binNum++){
            if(bins[binNum] != null){
                bins[binNum].calcFrequency(total);
            }
        }
    }

   public void addToBin(int binNum){
       bins[binNum].add();
       numTotalValues++;
   }

    // private functions

    private void createBins(int numBins){
        bins = new Bin[numBins];
        int     binId        = 0;
        double  range        = max - min;
        double  lower        = min;

        binInterval          = range / numBins;

        while(binId < numBins){
            bins[binId] = new Bin(lower, lower+binInterval);
            lower += binInterval;
            binId++;
        }
    }


}
