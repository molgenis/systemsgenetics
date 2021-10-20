package betaqtl.stat;

import JSci.maths.ArrayMath;
import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.tint.IntComparator;
import org.apache.commons.collections.primitives.ArrayDoubleList;
import org.apache.commons.math3.stat.ranking.NaNStrategy;
import org.apache.commons.math3.stat.ranking.NaturalRanking;
import org.apache.commons.math3.stat.ranking.RankingAlgorithm;
import org.apache.commons.math3.stat.ranking.TiesStrategy;

import java.util.Arrays;
import java.util.HashSet;

public class RankArray {
    private static final RankingAlgorithm COV_RANKER_TIE;
    private static final RankingAlgorithm COV_RANKER;
    public double[] xdouble = null;
    public int[] ydouble = null;
    public Swapper swapperdouble = null;
    public IntComparator compdouble = null;
    public float[] x = null;
    public int[] y = null;
    public Swapper swapper = null;
    public IntComparator comp = null;

    public RankArray() {
        this.swapperdouble = new Swapper() {
            public void swap(int a, int b) {
                double t1 = RankArray.this.xdouble[a];
                RankArray.this.xdouble[a] = RankArray.this.xdouble[b];
                RankArray.this.xdouble[b] = t1;
                int t2 = RankArray.this.ydouble[a];
                RankArray.this.ydouble[a] = RankArray.this.ydouble[b];
                RankArray.this.ydouble[b] = t2;
            }
        };
        this.compdouble = new IntComparator() {
            public int compare(int a, int b) {
                return RankArray.this.xdouble[a] == RankArray.this.xdouble[b] ? 0 : (RankArray.this.xdouble[a] < RankArray.this.xdouble[b] ? -1 : 1);
            }
        };
        this.swapper = new Swapper() {
            public void swap(int a, int b) {
                float t1 = RankArray.this.x[a];
                RankArray.this.x[a] = RankArray.this.x[b];
                RankArray.this.x[b] = t1;
                int t2 = RankArray.this.y[a];
                RankArray.this.y[a] = RankArray.this.y[b];
                RankArray.this.y[b] = t2;
            }
        };
        this.comp = new IntComparator() {
            public int compare(int a, int b) {
                return RankArray.this.x[a] == RankArray.this.x[b] ? 0 : (RankArray.this.x[a] < RankArray.this.x[b] ? -1 : 1);
            }
        };
    }

    public double[] rank(double[] x, boolean giveTiesSameRank) {
        double[] rankNaN = new double[x.length];
        // check if input is all NaN
        int nanCtr = 0;
        for (int v = 0; v < x.length; v++) {
            if (Double.isNaN(x[v])) {
                nanCtr++;
            }
        }
        if (nanCtr == x.length) {
            Arrays.fill(rankNaN, Double.NaN);
            return rankNaN;
        }

        double[] rank;
        if (!giveTiesSameRank) {
            rank = COV_RANKER.rank(x);
        } else {
            rank = COV_RANKER_TIE.rank(x);
        }


        int vctr = 0;
        for (int v = 0; v < x.length; ++v) {
            if (Double.isNaN(x[v])) {
                rankNaN[v] = Double.NaN;
            } else {
                rankNaN[v] = rank[vctr] - 1;
                vctr++;
            }

        }

        return rankNaN;
    }

    public float[] rank(float[] x, boolean giveTiesSameRank) {
        this.x = (float[]) x.clone();
        this.y = new int[x.length];

        for (int v = 0; v < x.length; this.y[v] = v++) {
        }

        GenericSorting.quickSort(0, x.length, this.comp, this.swapper);
        float[] rank = new float[x.length];

        for (int v = 0; v < x.length; ++v) {
            rank[this.y[v]] = (float) v;
        }

        if (!giveTiesSameRank) {
            return rank;
        } else {
            this.fixTiesFloat(rank, x);
            return rank;
        }
    }

    private void fixTiesFloat(float[] rank, float[] x) {
        HashSet<Float> fixedValues = new HashSet();

        for (int i = 0; i < x.length; ++i) {
            for (int j = i + 1; j < x.length; ++j) {
                if (x[i] == x[j] && !fixedValues.contains(x[i])) {
                    this.replaceRankFloat(x[i], i, j, rank, x);
                    fixedValues.add(x[i]);
                    break;
                }
            }
        }

    }

    private void replaceRankFloat(float f, int i, int j, float[] rank, float[] x) {
        ArrayDoubleList t = new ArrayDoubleList();
        t.add((double) rank[i]);
        t.add((double) rank[j]);

        for (int k = j + 1; k < x.length; ++k) {
            if (x[k] == f) {
                t.add((double) rank[k]);
            }
        }

        double newRank = ArrayMath.mean(t.toArray());
        rank[i] = (float) newRank;
        rank[j] = (float) newRank;

        for (int k = j + 1; k < x.length; ++k) {
            if (x[k] == f) {
                rank[k] = (float) newRank;
            }
        }

    }

    static {
        COV_RANKER_TIE = new NaturalRanking(NaNStrategy.REMOVED, TiesStrategy.AVERAGE);
        COV_RANKER = new NaturalRanking(NaNStrategy.REMOVED, TiesStrategy.SEQUENTIAL);
    }
}
