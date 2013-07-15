/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats;

import JSci.maths.ArrayMath;
import JSci.maths.statistics.NormalDistribution;
import java.io.IOException;
import java.util.Arrays;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author juha
 */
public class FDRUtils {

    private static NormalDistribution normDist;

    public enum CloneArrays {

        CLONE_BOTH, DONT_CLONE, CLONE_REAL
    }

    public enum SortArrays {

        SORT_DESCENDING, ALREADY_SORTED_DESCENDING
    }

    public enum AbsoluteValues {

        USE_ABSOLUTE, USE_AS_IS
    }

    // prevent instantiation, only use static factory methods
    private FDRUtils() {
        normDist = null;
    }

    /**
     *
     * Returns the threshold value in the whole of given real data matrix that
     * corresponds to the given false discovery rate against the null
     * distribution from given permuted data matrix.
     *
     * @param real real data [x][y]
     * @param permuted permuted data [x][y]
     * @param fdr false discovery rate to find
     * @param c should the data be cloned or not (because the arrays are
     * modified)
     * @param s should the data be sorted first or are they already sorted
     * @param a should absolute values be used
     * @return largest threshold value in real data so that FDR is not larger
     * than the given one
     */
    public static double getGlobalThreshold(double[][] real, double[][] permuted, double fdr, CloneArrays c, SortArrays s, AbsoluteValues a) {

        double[][][] perm3d = new double[1][][];
        perm3d[0] = permuted;
        return getGlobalThreshold(real, perm3d, fdr, c, s, a);
    }

    /**
     *
     * Calculates the number of significant Z-scores per each row in the given Z-score matrix. Z-scores are converted to p-values and the given p-value threshold is used as cutoff for significance. A two-tailed test is assumed.
     *
     * @param data a Z-score matrix
     * @param threshold p-value threshold to use
     * @param s should the data be sorted first or are they already sorted
     * @param a should absolute values be used
     * @return an array with the number of significant (according to the given p-value threshold) values in the matrix per row
     */
    public static int[] getNrSignificant(double[][] data, double threshold, SortArrays s, AbsoluteValues a) {

        switch (s) {
            case SORT_DESCENDING:
                data = sortDesc(data);
                break;
        }

        switch (a) {
            case USE_ABSOLUTE:
                data = ArrayMath.abs(data);
                break;
        }

        int[] nrSignificant = new int[data.length];
        for (int gene = 0; gene < data.length; gene++) {
            for (int term = 0; term < data[0].length; term++) {
                if (zToP(data[gene][term]) < threshold) {
                    nrSignificant[gene]++;
                } else {
                    break;
                }
            }
        }

        return nrSignificant;
    }

    /**
     *
     * Returns the p-value corresponding to the given Z-score assuming a two tailed-test.
     *
     * @param z A Z-score to convert
     * @return The p-value corresponding to the Z-score
     */
    public static double zToP(double z) {

        if (normDist == null) {
            normDist = new NormalDistribution();
        }
        double p;
        if (z > 0) {
            p = normDist.cumulative(-z);
        } else {
            p = normDist.cumulative(z);
        }
        if (p > 0.5) {
            p = 1 - p;
        }
        p *= 2.0d;

        return p;
    }

    /**
     *
     * Returns the threshold value in the whole of given real data matrix that
     * corresponds to the given false discovery rate against the null
     * distribution from given permuted data matrices.
     *
     * @param real real data [x][y]
     * @param permuted permuted data [permutation][x][y]
     * @param fdr false discovery rate to find
     * @param c should the data be cloned or not (because the arrays are
     * modified)
     * @param s should the data be sorted first or are they already sorted
     * @param a should absolute values be used
     * @return largest threshold value in real data so that FDR is not larger
     * than the given one
     */
    public static double getGlobalThreshold(double[][] real, double[][][] permuted, double fdr, CloneArrays c, SortArrays s, AbsoluteValues a) {

        checkNulls(real, permuted);
        checkSizes(real, permuted);
        checkFDR(fdr);

        int nrPermutations = permuted.length;
        double[][] realData = real;
        double[][][] permData = permuted;
        switch (c) {
            case CLONE_BOTH:
                realData = real.clone();
                permData = permuted.clone();
                break;
            case CLONE_REAL:
                realData = real.clone();
                break;
        }

        switch (s) {
            case SORT_DESCENDING:
                realData = sortDesc(realData);
                permData = sortDesc(permData);
                break;
        }

        switch (a) {
            case USE_ABSOLUTE:
                realData = ArrayMath.abs(realData);
                for (int p = 0; p < nrPermutations; p++) {
                    permData[p] = ArrayMath.abs(permData[p]);
                }
                break;
        }

        double threshold = -1;
        int nrSignificantInReal = 0;
        int nrSignificantInPerm = 0;
        int[] realColsTraversed = new int[realData.length];
        int[][] permColsTraversed = new int[nrPermutations][realData.length];
        double foundFDR = -1;

        while (foundFDR < fdr) {

            // get next largest value i.e. threshold from real data
            // knowing that the data are sorted both by rows and columns
            double realThreshold = 0;
            int nextMaxRow = -1;
            for (int i = 0; i < realColsTraversed.length; i++) {
                if (realData[i][realColsTraversed[i]] > realThreshold) {
                    realThreshold = realData[i][realColsTraversed[i]];
                    nextMaxRow = i;
                }
                if (realColsTraversed[i] == 0) {
                    break;
                }
            }
            realColsTraversed[nextMaxRow]++;
            nrSignificantInReal++;

            for (int p = 0; p < nrPermutations; p++) {
                for (int i = 0; i < permColsTraversed[0].length; i++) {
                    while (permData[p][i][permColsTraversed[p][i]] > realThreshold) {
                        nrSignificantInPerm++;
                        permColsTraversed[p][i]++;
                    }
                    if (permColsTraversed[p][i] == 0) {
                        break;
                    }
                }
            }

            foundFDR = (double) nrSignificantInPerm / (double) nrPermutations / (double) nrSignificantInReal;
            if (foundFDR <= fdr) {
                threshold = realThreshold;
            }
        }

        return threshold;
    }

    /**
     *
     * Returns the threshold values per row in the given real data matrix that
     * corresponds to the given false discovery rate against the null
     * distribution from rows of given permuted data matrix.
     *
     * @param real real data [x][y]
     * @param permuted permuted data [x][y]
     * @param fdr false discovery rate to find
     * @param c should the data be cloned or not (because the arrays are
     * modified)
     * @param s should the data be sorted first or are they already sorted
     * @param a should absolute values be used
     * @return largest threshold values per row in real data so that FDR is not
     * larger than the given one
     */
    public static double[] getThresholdsPerRow(double[][] real, double[][] permuted, double fdr, CloneArrays c, SortArrays s, AbsoluteValues a) {

        double[][][] perm3d = new double[1][][];
        perm3d[0] = permuted;
        return getThresholdsPerRow(real, perm3d, fdr, c, s, a);
    }

    /**
     *
     * Returns the threshold values per row in the given real data matrix that
     * corresponds to the given false discovery rate against the null
     * distribution from rows of given permuted data matrices.
     *
     * @param real real data [x][y]
     * @param permuted permuted data [permutation][x][y]
     * @param fdr false discovery rate to find
     * @param c should the data be cloned or not (because the arrays are
     * modified)
     * @param s should the data be sorted first or are they already sorted
     * @param a should absolute values be used
     * @return largest threshold values per row in real data so that FDR is not
     * larger than the given one
     */
    public static double[] getThresholdsPerRow(double[][] real, double[][][] permuted, double fdr, CloneArrays c, SortArrays s, AbsoluteValues a) {

        checkNulls(real, permuted);
        checkSizes(real, permuted);
        checkFDR(fdr);

        int nrPermutations = permuted.length;
        double[][] realData = real;
        double[][][] permData = permuted;
        switch (c) {
            case CLONE_BOTH:
                realData = real.clone();
                permData = permuted.clone();
                break;
            case CLONE_REAL:
                realData = real.clone();
                break;
        }

        switch (a) {
            case USE_ABSOLUTE:
                realData = ArrayMath.abs(realData);
                for (int p = 0; p < nrPermutations; p++) {
                    permData[p] = ArrayMath.abs(permData[p]);
                }
                break;
        }

        double[] thresholds = new double[realData.length];

        for (int row = 0; row < realData.length; row++) {

            if (s == SortArrays.SORT_DESCENDING) {
                Arrays.sort(realData[row]);
                realData[row] = ArrayMath.invert(realData[row]);
                for (int p = 0; p < nrPermutations; p++) {
                    Arrays.sort(permData[p][row]);
                    permData[p][row] = ArrayMath.invert(permData[p][row]);
                }
            }

            double[] thresholdPerm = new double[nrPermutations];
            double foundFDR = 1;

            for (int p = 0; p < nrPermutations; p++) {
                int nrSignificantInReal = 0;
                thresholdPerm[p] = realData[row][0];
                for (int i = 0; i < realData[0].length; i++) {
                    double currentThreshold = realData[row][i];
                    nrSignificantInReal++;
                    int nrSignificantInPerm = 0;
                    for (int j = 0; j < permData[0][0].length; j++) {
                        if (permData[p][row][j] > currentThreshold) {
                            nrSignificantInPerm++;
                        } else {
                            break;
                        }
                    }

                    foundFDR = nrSignificantInPerm / nrSignificantInReal;

                    if (foundFDR > fdr) {
                        break;
                    } else {
                        thresholdPerm[p] = currentThreshold;
                    }
                }
            }

            Arrays.sort(thresholdPerm);
            thresholds[row] = thresholdPerm[(int) (0.8 * (double) nrPermutations) - 1];

        }

        return thresholds;
    }

    public static int[] getNrSignificant(double[][] data, double[] thresholds, CloneArrays c, SortArrays s, AbsoluteValues a) {

        int[] nrSignificant = new int[thresholds.length];
        double[][] realData = data;
        switch (c) {
            case CLONE_BOTH:
                realData = data.clone();
                break;
            case CLONE_REAL:
                realData = data.clone();
                break;
        }

        switch (a) {
            case USE_ABSOLUTE:
                realData = ArrayMath.abs(realData);
                break;
        }

        if (s == SortArrays.SORT_DESCENDING) {
            for (int row = 0; row < realData.length; row++) {
                Arrays.sort(realData[row]);
                realData[row] = ArrayMath.invert(realData[row]);
            }
        }

        for (int row = 0; row < realData.length; row++) {
            int nrSignificantThisRow = 0;
            for (int i = 0; i < realData[0].length; i++) {
                if (realData[row][i] >= thresholds[row]) {
                    nrSignificantThisRow++;
                }
            }
            nrSignificant[row] = nrSignificantThisRow;
        }
        return nrSignificant;
    }

    public static double[] getThresholdsPerRow(double[][] real, String[] permutedFiles, double fdr, double mvp, CloneArrays c, SortArrays s, AbsoluteValues a) throws IOException {

        checkFDR(fdr);

        int nrPermutations = permutedFiles.length;
        double[][] realData = real;
        switch (c) {
            case CLONE_BOTH:
                realData = real.clone();
                break;
            case CLONE_REAL:
                realData = real.clone();
                break;
        }

        switch (a) {
            case USE_ABSOLUTE:
                realData = ArrayMath.abs(realData);
                break;
        }

        double[] thresholds = new double[realData.length];

        if (s == SortArrays.SORT_DESCENDING) {
            for (int row = 0; row < realData.length; row++) {
                Arrays.sort(realData[row]);
                realData[row] = ArrayMath.invert(realData[row]);
            }
        }

        double[][] thresholdRowPerm = new double[realData.length][nrPermutations];

        for (int p = 0; p < nrPermutations; p++) {

            System.out.println("Processing permutation " + p);

            DoubleMatrixDataset<String, String> permDataset = new DoubleMatrixDataset<String, String>(permutedFiles[p]);
            double[][] permData = permDataset.getRawDataTransposed();
            if (a == AbsoluteValues.USE_ABSOLUTE) {
                permData = ArrayMath.abs(permData);
            }

            double foundFDR = 1;
            for (int row = 0; row < realData.length; row++) {

                if (s == SortArrays.SORT_DESCENDING) {
                    Arrays.sort(permData[row]);
                    permData[row] = ArrayMath.invert(permData[row]);
                }

                int nrSignificantInReal = 0;
                thresholdRowPerm[row][p] = Double.MAX_VALUE;
                for (int i = 0; i < realData[0].length; i++) {
                    double currentThreshold = realData[row][i];
                    nrSignificantInReal++;
                    int nrSignificantInPerm = 0;
                    for (int j = 0; j < permData[0].length; j++) {
                        if (permData[row][j] > currentThreshold) {
                            nrSignificantInPerm++;
                        } else {
                            break;
                        }
                    }
                    foundFDR = (double) nrSignificantInPerm / nrSignificantInReal;
//                    System.out.println(nrSignificantInReal + "\t" + nrSignificantInPerm + "\t" + foundFDR);
                    if (row == 16846) {
//                        System.out.println("ENSG00000129170 " + nrSignificantInReal + " " + nrSignificantInPerm + " " + foundFDR);
                    }
                    if (foundFDR > fdr) {
                        break;
                    } else {
                        thresholdRowPerm[row][p] = currentThreshold;
                    }
                }
            }
        }

        // 80 % confidence
        for (int row = 0; row < realData.length; row++) {
            Arrays.sort(thresholdRowPerm[row]);
            thresholds[row] = thresholdRowPerm[row][(int) (mvp * (double) nrPermutations) - 1];
        }

        return thresholds;
    }

    /**
     *
     * Sorts the given array both by rows and by columns.
     *
     * @param array array to be sorted
     * @return sorted array
     */
    private static double[][] sortDesc(double[][] array) {
        // sort data wrt each row
        for (int i = 0; i < array.length; i++) {
            Arrays.sort(array[i]);
            array[i] = ArrayMath.invert(array[i]);
        }
        // sort data wrt each column
        array = ArrayMath.transpose(array);
        for (int i = 0; i < array.length; i++) {
            Arrays.sort(array[i]);
            array[i] = ArrayMath.invert(array[i]);
        }
        array = ArrayMath.transpose(array);

//        for (int i = 0; i < 10; i++) {
//            for (int j = 0; j < 10; j++) {
//                System.out.print(array[i][j] + "\t");
//            }
//            System.out.println();
//        }

        return array;
    }

    /**
     *
     * Sorts the given arrays both by rows and by columns.
     *
     * @param array arrays to be sorted ([array][row][col])
     * @return sorted arrays
     */
    private static double[][][] sortDesc(double[][][] array) {

        for (int p = 0; p < array.length; p++) {
            // sort data wrt each row
            for (int i = 0; i < array[p].length; i++) {
                Arrays.sort(array[p][i]);
                array[p][i] = ArrayMath.invert(array[p][i]);
            }
            // sort data wrt each column
            array[p] = ArrayMath.transpose(array[p]);
            for (int i = 0; i < array[p].length; i++) {
                Arrays.sort(array[p][i]);
                array[p][i] = ArrayMath.invert(array[p][i]);
            }
            array[p] = ArrayMath.transpose(array[p]);
        }

        return array;
    }

    private static void checkNulls(double[][] real, double[][][] permuted) {
        if (real == null) {
            throw new IllegalArgumentException("Null real data given.");
        }
        if (permuted == null) {
            throw new IllegalArgumentException("Null permuted data given.");
        }
        for (int i = 0; i < permuted.length; i++) {
            if (permuted[i] == null) {
                throw new IllegalArgumentException("Null permuted data given at index " + i + ".");
            }
        }

    }

    private static void checkSizes(double[][] real, double[][][] permuted) {
        if (real.length != permuted[0].length || real[0].length != permuted[0][0].length) {
            throw new IllegalArgumentException("Matrix sizes (real and permuted) must be the same. Real: " + real.length + " x " + real[0].length
                    + ", permuted: " + permuted[0].length + " x " + permuted[0][0].length);
        }
    }

    private static void checkFDR(double fdr) {
        if (fdr < Double.MIN_VALUE || fdr >= 1) {
            throw new IllegalArgumentException("FDR must be larger than 0 and smaller than 1: " + fdr);
        }
    }
}
