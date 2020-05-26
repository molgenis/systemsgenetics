/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.matrix2;

import cern.colt.function.tint.IntComparator;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;

/**
 * This is identical to code in cern.colt.matrix.tdouble.algo.DoubleSorting but now public static
 * 
 * 
 * @author patri
 */
public class DoubleMatrix1dOrder {

	/**
	 * Sorts indexes of the <code>vector</code> into ascending order.
	 *
	 * @param vector
	 * @return sorted indexes
	 */
	public static int[] sortIndex(final DoubleMatrix1D vector) {
		int[] indexes = new int[(int) vector.size()]; // row indexes to reorder
		// instead of matrix itself
		for (int i = indexes.length; --i >= 0;) {
			indexes[i] = i;
		}
		IntComparator comp = null;
		if (vector instanceof DenseDoubleMatrix1D) {
			final double[] velems = (double[]) vector.elements();
			final int zero = (int) vector.index(0);
			final int stride = vector.stride();
			comp = new IntComparator() {
				public int compare(int a, int b) {
					int idxa = zero + a * stride;
					int idxb = zero + b * stride;
					double av = velems[idxa];
					double bv = velems[idxb];
					if (av != av || bv != bv) {
						return compareNaN(av, bv); // swap NaNs to the end
					}
					return av < bv ? -1 : (av == bv ? 0 : 1);
				}
			};
		} else {
			comp = new IntComparator() {
				public int compare(int a, int b) {
					double av = vector.getQuick(a);
					double bv = vector.getQuick(b);
					if (av != av || bv != bv) {
						return compareNaN(av, bv); // swap NaNs to the end
					}
					return av < bv ? -1 : (av == bv ? 0 : 1);
				}
			};
		}

		runSort(indexes, 0, indexes.length, comp);

		return indexes;
	}
	
	/**
	 * Sorts indexes of the <code>vector</code> into decending order.
	 * 
	 * Adapted form sortIndex
	 *
	 * @param vector
	 * @return sorted indexes
	 */
	public static int[] sortIndexReverse(final DoubleMatrix1D vector) {
		int[] indexes = new int[(int) vector.size()]; // row indexes to reorder
		// instead of matrix itself
		for (int i = indexes.length; --i >= 0;) {
			indexes[i] = i;
		}
		IntComparator comp = null;
		if (vector instanceof DenseDoubleMatrix1D) {
			final double[] velems = (double[]) vector.elements();
			final int zero = (int) vector.index(0);
			final int stride = vector.stride();
			comp = new IntComparator() {
				public int compare(int a, int b) {
					int idxa = zero + a * stride;
					int idxb = zero + b * stride;
					double av = velems[idxa];
					double bv = velems[idxb];
					if (av != av || bv != bv) {
						return compareNaN(av, bv); // swap NaNs to the end
					}
					return av > bv ? -1 : (av == bv ? 0 : 1);
				}
			};
		} else {
			comp = new IntComparator() {
				public int compare(int a, int b) {
					double av = vector.getQuick(a);
					double bv = vector.getQuick(b);
					if (av != av || bv != bv) {
						return compareNaN(av, bv); // swap NaNs to the end
					}
					return av > bv ? -1 : (av == bv ? 0 : 1);
				}
			};
		}

		runSort(indexes, 0, indexes.length, comp);

		return indexes;
	}

	/**
	 * Compare two values, one of which is assumed to be Double.NaN
	 */
	private static int compareNaN(double a, double b) {
		if (a != a) {
			if (b != b) {
				return 0; // NaN equals NaN
			} else {
				return 1; // e.g. NaN > 5
			}
		}
		return -1; // e.g. 5 < NaN
	}
	
	protected static void runSort(int[] a, int fromIndex, int toIndex, IntComparator c) {
        cern.colt.Sorting.parallelQuickSort(a, fromIndex, toIndex, c);
    }

}
