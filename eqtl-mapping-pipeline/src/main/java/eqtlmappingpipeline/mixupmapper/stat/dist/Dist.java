/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package eqtlmappingpipeline.mixupmapper.stat.dist;

import java.util.Collections;
import java.util.ArrayList;
/**
 *
 * @author harm-jan
 */
public class Dist {
    protected double max, min;
    protected double var, se, mean, median;
  
    public Dist(){
        max     = 0.0;
        min     = 0.0;
        mean    = 0.0;
        median  = 0.0;
        var     = 0.0;
        se      = 0.0;  
    }

    public void createDist(ArrayList<Double> alValues){
        max = Collections.max(alValues);
        min = Collections.min(alValues);

        calcStatistics(alValues);
    }

    public double getMax(){
        return max;
    }

    public double getMin(){
        return min;
    }

    public double getVariance(){
        return var;
    }

    public void setVariance(double newVar){
        var = newVar;
    }

    public double getSD(){
        return Math.sqrt(var);
    }

    public double getMedian(){
        return median;
    }

    public void setMedian(double newMedian){
        median = newMedian;
    }

    public double getMean(){
        return mean;
    }

    public void setMean(double newMean){
        mean = newMean;
    }

// PROTECTED functions
    protected void calcStatistics(ArrayList<Double> alValues){

        double total = 0.0;
        for(Double val: alValues){
           total += val;
        }

        int size = alValues.size();

        mean        = total / size;
        var         = 0.0;

        for(Object val: alValues){
           var     += ((Double) val - mean) * ((Double) val - mean);
        }

        Collections.sort(alValues);

        median = 0.0;
        if(size % 2 == 0) {
           // even number
            int mediantop   = (int) (size / 2 + 0.5);
            int medianbottom= (int) (size / 2 - 0.5);
            median	    = ((Double) alValues.get(mediantop) + (Double) alValues.get(medianbottom) )/ 2;
        } else {
            int index	    = (int) (size / 2 + 0.5);
            median	    = (Double) alValues.get(index);
        }
    }
    
    
}
