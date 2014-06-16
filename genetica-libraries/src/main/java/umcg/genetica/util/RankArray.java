/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.util;

import cern.colt.GenericSorting;
import java.util.HashSet;
import org.apache.commons.collections.primitives.ArrayDoubleList;
import cern.colt.function.tint.IntComparator;
import cern.colt.Swapper;

/**
 *
 * @author harmjan
 */
public class RankArray {
//    public static double[] rank(double[] x){
//        umcg.genetica.util.Rank rank = new umcg.genetica.util.Rank(x, 0d);
//        double[] ranks =  rank.getRanks();
//        double[] returnRanks = new double[ranks.length];
//        for(int i=0; i<returnRanks.length; i++){
//            returnRanks[i] = ranks[i] - 1;
//        }
//        return returnRanks;
//    }
//    
//    public static float[] rank(float[] x){
//        double[] data = new double[x.length];
//        for(int i=0; i<data.length; i++){
//            data[i] = (double) x[i];
//        }
//        jsc.util.Rank rank = new jsc.util.Rank(data, 0d);
//        double[] ranks = rank.getRanks();
//        float[] dranks = new float[ranks.length];
//        for(int i=0; i<dranks.length; i++){
//            dranks[i] = (float) ranks[i] - 1;
//        }
//        return dranks;
//    }

    public double[] xdouble = null;
    public int[] ydouble = null;
    public Swapper swapperdouble = null;
    public IntComparator compdouble = null;
    public float[] x = null;
    public int[] y = null;
    public Swapper swapper = null;
    public IntComparator comp = null;

    public RankArray() {
        swapperdouble = new cern.colt.Swapper() {
            @Override
            public void swap(int a, int b) {
                double t1;
                int t2;
                t1 = xdouble[a];
                xdouble[a] = xdouble[b];
                xdouble[b] = t1;
                t2 = ydouble[a];
                ydouble[a] = ydouble[b];
                ydouble[b] = t2;
            }
        };

        compdouble = new IntComparator() {
            @Override
            public int compare(int a, int b) {
                return xdouble[a] == xdouble[b] ? 0 : (xdouble[a] < xdouble[b] ? -1 : 1);
            }
        };

        swapper = new cern.colt.Swapper() {
            @Override
            public void swap(int a, int b) {
                float t1;
                int t2;
                t1 = x[a];
                x[a] = x[b];
                x[b] = t1;
                t2 = y[a];
                y[a] = y[b];
                y[b] = t2;
            }
        };

        comp = new IntComparator() {
            @Override
            public int compare(int a, int b) {
                return x[a] == x[b] ? 0 : (x[a] < x[b] ? -1 : 1);
            }
        };
    }

    public double[] rank(double[] x, boolean giveTiesSameRank) {
        
        this.xdouble = x.clone();
        ydouble = new int[x.length];
        for (int v = 0; v < x.length; v++) {
            ydouble[v] = v;
        }
        GenericSorting.quickSort(0, x.length, compdouble, swapperdouble);
        double[] rank = new double[x.length];
        for (int v = 0; v < x.length; v++) {
            rank[ydouble[v]] = v;
        }
        
        if (!giveTiesSameRank) {
            return rank;
        } else {
            fixTiesDouble(rank, x);
            return rank;
        }
    }

    public float[] rank(float[] x, boolean giveTiesSameRank) {
        
        this.x = x.clone();
        y = new int[x.length];
        for (int v = 0; v < x.length; v++) {
            y[v] = v;
        }
        GenericSorting.quickSort(0, x.length, comp, swapper);
        float[] rank = new float[x.length];
        for (int v = 0; v < x.length; v++) {
            rank[y[v]] = v;
        }
        if (!giveTiesSameRank) {
            return rank;
        } else {
            fixTiesFloat(rank, x);
            return rank;
        }
    }

    private void fixTiesDouble(double[] rank, double[] x) {
        HashSet<Double> fixedValues = new HashSet<Double>();
        
        for(int i=0; i<x.length;++i){
            for(int j=i+1; j<x.length;++j){
                if(x[i] == x[j] && !fixedValues.contains(x[i])){
                    replaceRankDouble(x[i], i, j, rank, x);
                    fixedValues.add(x[i]);
                    break;
                }
            }
        }
        
    }
    
    private void replaceRankDouble(double d, int i, int j, double[] rank, double[] x) {
        ArrayDoubleList t = new ArrayDoubleList();
        
        t.add(rank[i]);
        t.add(rank[j]);
        
        for(int k = j+1; k < x.length; k++){
            if(x[k]==d){
                t.add(rank[k]);
            }
        }
        double newRank = JSci.maths.ArrayMath.mean(t.toArray());
        
        rank[i] = newRank;
        rank[j] = newRank;
        
        for(int k = j+1; k < x.length; k++){
            if(x[k]==d){
                rank[k] = newRank;
            }
        }
    }
    
    private void fixTiesFloat(float[] rank, float[] x) {
        HashSet<Float> fixedValues = new HashSet<Float>();
        
        for(int i=0; i<x.length;++i){
            for(int j=i+1; j<x.length;++j){
                if(x[i] == x[j] && !fixedValues.contains(x[i])){
                    replaceRankFloat(x[i], i, j, rank, x);
                    fixedValues.add(x[i]);
                    break;
                }
            }
        }
    }

    private void replaceRankFloat(float f, int i, int j, float[] rank, float[] x) {
        ArrayDoubleList t = new ArrayDoubleList();
        
        t.add(rank[i]);
        t.add(rank[j]);
        
        for(int k = j+1; k < x.length; k++){
            if(x[k]==f){
                t.add(rank[k]);
            }
        }
        
        double newRank = JSci.maths.ArrayMath.mean(t.toArray());
        
        rank[i] = (float) newRank;
        rank[j] = (float) newRank;
        
        for(int k = j+1; k < x.length; k++){
            if(x[k]==f){
                rank[k] = (float) newRank;
            }
        }
    }
}
