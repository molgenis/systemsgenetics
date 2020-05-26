/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.util;

import cern.colt.GenericSorting;
import java.util.HashSet;
import org.apache.commons.collections.primitives.ArrayDoubleList;
import org.apache.commons.math3.stat.ranking.RankingAlgorithm;
import cern.colt.function.tint.IntComparator;
import cern.colt.Swapper;
import org.apache.commons.math3.stat.ranking.NaNStrategy;
import org.apache.commons.math3.stat.ranking.NaturalRanking;
import org.apache.commons.math3.stat.ranking.TiesStrategy;

/**
 *
 * @author harmjan
 */

//Please use the math3 natural ranker, especialy in tie resolving it is much faster!
//Note we can also chose to make this depend on natural ranker and do some more speed tweaks, no NaN strategy and only one tie fixer.
//Also make it directly 0 based in this case!

@Deprecated
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
    private static final RankingAlgorithm COV_RANKER_TIE = new NaturalRanking(NaNStrategy.FAILED, TiesStrategy.AVERAGE);
    private static final RankingAlgorithm COV_RANKER = new NaturalRanking(NaNStrategy.FAILED, TiesStrategy.SEQUENTIAL);
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
        double[] rank = null;
        if (!giveTiesSameRank) {
            rank = COV_RANKER.rank(x);
        } else {
            rank = COV_RANKER_TIE.rank(x);
        }
        for (int v = 0; v < rank.length; v++) {
                rank[v] = rank[v]-1;
            }
            return rank;
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

//    private void fixTiesDouble(double[] rank, double[] x) {
//        HashSet<Double> fixedValues = new HashSet<Double>();
//        
//        for(int i=0; i<x.length;++i){
//            for(int j=i+1; j<x.length;++j){
//                if(x[i] == x[j] && !fixedValues.contains(x[i])){
//                    replaceRankDouble(x[i], i, j, rank, x);
//                    fixedValues.add(x[i]);
//                    break;
//                }
//            }
//        }
//        
//    }
//    
//    private void replaceRankDouble(double d, int i, int j, double[] rank, double[] x) {
//        ArrayDoubleList t = new ArrayDoubleList();
//        
//        t.add(rank[i]);
//        t.add(rank[j]);
//        
//        for(int k = j+1; k < x.length; k++){
//            if(x[k]==d){
//                t.add(rank[k]);
//            }
//        }
//        double newRank = JSci.maths.ArrayMath.mean(t.toArray());
//        
//        rank[i] = newRank;
//        rank[j] = newRank;
//        
//        for(int k = j+1; k < x.length; k++){
//            if(x[k]==d){
//                rank[k] = newRank;
//            }
//        }
//    }
    
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
