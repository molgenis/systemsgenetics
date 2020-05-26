/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package eqtlmappingpipeline.mixupmapper.stat.dist;

/**
 *
 * @author harm-jan
 */
public class Bin {

    protected double    lower, upper;
    protected int       count;
    protected Double    frequency;
    protected int       maxCount;
    protected int       cumulative;
    protected Double    cumulativeFrequency;


    public Bin(double _lower, double _upper){
        lower               = _lower;
        upper               = _upper;
        frequency           = null;
        count               = 0;
        maxCount            = 0;
        cumulative          = 0;
        cumulativeFrequency = null;
    }

    public double getLowerBoundary(){
        return lower;
    }

    public double getUpperBoundary(){
        return upper;
    }

    public void setMaxCount(int max){
        maxCount = max;
     }

     public int getMaxCount(){
         return maxCount;
     }

     public void add(){
         count++;
     }

     public int getCount(){
         return count;
     }

     public int getCumulative(){
        return cumulative;
    }

    public void setCumulative(int sum){
        cumulative = sum;
    }

    public void calcFrequency(int total){
        frequency           = (double) count        / total;
        cumulativeFrequency = (double) cumulative   / total;
    }

    public double getFrequency(){
        return frequency;
    }

    public double getCumulativeFrequency(){
        return cumulativeFrequency;
    }
}
