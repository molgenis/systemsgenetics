/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.util;

/**
 * As taken from jsc.jar
 *
 * @author MarcJan
 */
public class Rank {

    private int n;
    private int s;
    private int t;
    private final int[] r;
    private final double[] rank;

    public Rank(double[] paramArrayOfDouble, double paramDouble) {
        this(paramArrayOfDouble.length, paramArrayOfDouble, paramDouble);
    }

    public Rank(int paramInt, double[] paramArrayOfDouble, double paramDouble) {
        if (paramInt < 1) {
            throw new IllegalArgumentException("No data to rank");
        }
        this.n = paramInt;

        double[] arrayOfDouble1 = new double[paramInt];
        this.r = new int[paramInt];
        for (int j = 0; j < paramInt; ++j) {
            this.r[j] = j;
        }

        this.s = 0;
        this.t = 0;

        double[] arrayOfDouble2 = new double[paramInt];

        System.arraycopy(paramArrayOfDouble, 0, arrayOfDouble2, 0, paramInt);

        Sort.sort(arrayOfDouble2, this.r, 0, paramInt - 1, true);

        arrayOfDouble1[(paramInt - 1)] = (paramInt - 1);
        paramInt--;
        for (int j = 0; j < paramInt; ++j) {
            if (Math.abs(arrayOfDouble2[j] - arrayOfDouble2[(j + 1)]) > paramDouble) {
                arrayOfDouble1[j] = j;
            } else {
                int k = 1;
                for (int i = j + 1; i < paramInt; ++i) {
                    if (Math.abs(arrayOfDouble2[i] - arrayOfDouble2[(i + 1)]) > paramDouble) {
                        break;
                    }
                    k++;
                }

                double d = j + 0.5D * k;
                for (int i1 = 0; i1 <= k; ++i1) {
                    arrayOfDouble1[(j + i1)] = d;
                }
                int m = k * (k + 1);
                this.s += m;
                this.t += m * (k + 2);
                j += k;
            }
        }
        paramInt++;
        this.rank = new double[paramInt];

        for (int i = 0; i < paramInt; ++i) {
            this.rank[this.r[i]] = (arrayOfDouble1[i] + 1.0D);
        }
    }

    public int getCorrectionFactor1() {
        return this.t;
    }

    public int getCorrectionFactor2() {
        return this.s;
    }

    public int getN() {
        return this.n;
    }

    public double getRank(int paramInt) {
        return this.rank[paramInt];
    }

    public double[] getRanks() {
        return this.rank;
    }

    public int getSortIndex(int paramInt) {
        return this.r[paramInt];
    }

    public int[] getSortIndexes() {
        return this.r;
    }

    public boolean hasTies() {
        return this.t > 0;
    }

}

