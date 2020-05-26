/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.matrix;

import umcg.genetica.containers.Pair;

/**
 *
 * @author juha
 */
public class MatrixTools {

    public static boolean containsZeros(double[][] a) {
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[0].length; j++) {
                if (a[i][j] == 0) {
                    return true;
                }
            }
        }
        return false;
    }

    public static boolean containsNaNs(double[][] a) {
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[0].length; j++) {
                if (Double.isNaN(a[i][j])) {
                    return true;
                }
            }
        }
        return false;
    }

    public static Pair<Integer, Integer> getIndexOfFirstZero(double[][] a) {
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[0].length; j++) {
                if (a[i][j] == 0) {
                    return new Pair<Integer, Integer>(i, j);
                }
            }
        }
        return null;
    }

    public static Pair<Integer, Integer> getIndexOfFirstNegative(double[][] a) {
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[0].length; j++) {
                if (a[i][j] < 0) {
                    return new Pair<Integer, Integer>(i, j);
                }
            }
        }
        return null;
    }

    public static Pair<Integer, Integer> getIndexOfFirstZeroIgnoreDiagonal(double[][] a) {
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[0].length; j++) {
                if (i != j && a[i][j] == 0) {
                    return new Pair<Integer, Integer>(i, j);
                }
            }
        }
        return null;
    }

    public static Pair<Integer, Integer> getIndexOfFirstNaN(double[][] a) {
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[0].length; j++) {
                if (Double.isNaN(a[i][j])) {
                    return new Pair<Integer, Integer>(i, j);
                }
            }
        }
        return null;
    }
}
