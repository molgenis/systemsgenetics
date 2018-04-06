/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.matrix2;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import no.uib.cipr.matrix.DenseMatrix;
import org.apache.commons.math3.stat.ranking.NaNStrategy;
import org.apache.commons.math3.stat.ranking.NaturalRanking;
import org.apache.commons.math3.stat.ranking.RankingAlgorithm;
import org.apache.commons.math3.stat.ranking.TiesStrategy;
import umcg.genetica.containers.Pair;
import umcg.genetica.math.stats.Descriptives;

/**
 *
 * @author MarcJan
 */
public class MatrixTools {

    public static void centerAndScaleColum(DoubleMatrix2D mat) {
        System.out.println("Standardizing probe mean and standard deviation");
        for (int c = 0; c < mat.columns(); c++) {
            double[] t = mat.viewColumn(c).toArray();
            double mean = Descriptives.mean(t);
            double stdev = Math.sqrt(Descriptives.variance(t, mean));
            for (int r = 0; r < mat.rows(); r++) {
                mat.set(r, c, ((mat.getQuick(r, c) - mean) / stdev));
            }
        }
    }

    public static void centerColum(DoubleMatrix2D mat) {
        System.out.println("Standardizing probe mean");
        for (int c = 0; c < mat.columns(); c++) {
            double mean = Descriptives.mean(mat.viewColumn(c).toArray());
            for (int r = 0; r < mat.rows(); r++) {
                mat.set(r, c, ((mat.getQuick(r, c) - mean)));
            }
        }
    }

    public static void scaleColum(DoubleMatrix2D mat) {
        System.out.println("Standardizing probe standard deviation");
        for (int c = 0; c < mat.columns(); c++) {
            double[] t = mat.viewColumn(c).toArray();
            double mean = Descriptives.mean(t);
            double stdev = Math.sqrt(Descriptives.variance(t, mean));
            for (int r = 0; r < mat.rows(); r++) {
                mat.set(r, c, ((mat.getQuick(r, c)) / stdev));
            }
        }
    }

    public static boolean containsZeros(DoubleMatrix2D a) {
        for (int i = 0; i < a.rows(); i++) {
            for (int j = 0; j < a.columns(); j++) {
                if (a.get(i, j) == 0) {
                    return true;
                }
            }
        }
        return false;
    }

    public static boolean containsNaNs(DoubleMatrix2D a) {
        for (int i = 0; i < a.rows(); i++) {
            for (int j = 0; j < a.columns(); j++) {
                if (Double.isNaN(a.get(i, j))) {
                    return true;
                }
            }
        }
        return false;
    }

    public static Pair<Integer, Integer> getIndexOfFirstZero(DoubleMatrix2D a) {
        for (int i = 0; i < a.rows(); i++) {
            for (int j = 0; j < a.columns(); j++) {
                if (a.get(i, j) == 0) {
                    return new Pair<Integer, Integer>(i, j);
                }
            }
        }
        return null;
    }

    public static Pair<Integer, Integer> getIndexOfFirstNegative(DoubleMatrix2D a) {
        for (int i = 0; i < a.rows(); i++) {
            for (int j = 0; j < a.columns(); j++) {
                if (a.get(i, j) < 0) {
                    return new Pair<Integer, Integer>(i, j);
                }
            }
        }
        return null;
    }

    public static Pair<Integer, Integer> getIndexOfFirstZeroIgnoreDiagonal(DoubleMatrix2D a) {
        for (int i = 0; i < a.rows(); i++) {
            for (int j = 0; j < a.columns(); j++) {
                if (i != j && a.get(i, j) == 0) {
                    return new Pair<Integer, Integer>(i, j);
                }
            }
        }
        return null;
    }

    public static Pair<Integer, Integer> getIndexOfFirstNaN(DoubleMatrix2D a) {
        for (int i = 0; i < a.rows(); i++) {
            for (int j = 0; j < a.columns(); j++) {
                if (Double.isNaN(a.get(i, j))) {
                    return new Pair<Integer, Integer>(i, j);
                }
            }
        }
        return null;
    }

    public static DenseDoubleMatrix2D toDenseDoubleMatrix(DenseMatrix matrix) {
        DenseDoubleMatrix2D matrix2D = new DenseDoubleMatrix2D(matrix.numRows(), matrix.numColumns());
        
        for(int c = 0; c<matrix.numColumns(); c++){
            for(int r = 0; r<matrix.numRows(); r++){
                matrix2D.setQuick(r,c, matrix.get(r, c));
            }
        }
        
        return matrix2D;
    }

    public static void rankColumns(DoubleMatrix2D matrix) {
        RankingAlgorithm COV_RANKER_TIE = new NaturalRanking(NaNStrategy.FAILED, TiesStrategy.AVERAGE);
        
        for (int c = 0; c < matrix.columns(); c++) {
            double[] rank = COV_RANKER_TIE.rank(matrix.viewColumn(c).toArray());
            for (int r = 0; r < matrix.rows(); r++) {
                matrix.set(r, c, (rank[r]));
            }
        }
    }
}
