/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package eqtlmappingpipeline.util;

/**
 *
 * @author harm-jan
 */
public class StringArrayByDoubleArraySort {

    public static void sort(double[] a, String[] b) {
        if(a.length == b.length){
            sort2(a, b, 0, a.length);
        }
    }

    private static void sort2(double a[], String[] b, int fromIndex, int toIndex) {
        final long NEG_ZERO_BITS = Double.doubleToLongBits(-0.0d);
        /*
         * The sort is done in three phases to avoid the expense of using
         * NaN and -0.0 aware comparisons during the main sort.
         */

        /*
         * Preprocessing phase:  Move any NaN's to end of array, count the
         * number of -0.0's, and turn them into 0.0's.
         */
        int numNegZeros = 0;
        int i = fromIndex, n = toIndex;
        while(i < n) {
            if (a[i] != a[i]) {
		double swap = a[i];
                String swapDesc = b[i];
                b[i] = b[--n];
                a[i] = a[--n];
                a[n] = swap;
                b[n] = swapDesc;
            } else {
                if (a[i]==0 && Double.doubleToLongBits(a[i])==NEG_ZERO_BITS) {
                    a[i] = 0.0d;
                    numNegZeros++;
                }
                i++;
            }
        }

        // Main sort phase: quicksort everything but the NaN's
	sort1(a, b, fromIndex, n-fromIndex);

        // Postprocessing phase: change 0.0's to -0.0's as required
        if (numNegZeros != 0) {
            int j = binarySearch0(a, fromIndex, n, 0.0d); // posn of ANY zero
            do {
                j--;
            } while (j>=0 && a[j]==0.0d);

            // j is now one less than the index of the FIRST zero
            for (int k=0; k<numNegZeros; k++)
                a[++j] = -0.0d;
        }
    }
    
    private static void sort1(double x[], String[] y, int off, int len) {
	// Insertion sort on smallest arrays
	if (len < 7) {
	    for (int i=off; i<len+off; i++)
		for (int j=i; j>off && x[j-1]>x[j]; j--) {
		    swap(x, y, j, j-1);
                }
	    return;
	}

	// Choose a partition element, v
	int m = off + (len >> 1);       // Small arrays, middle element
	if (len > 7) {
	    int l = off;
	    int n = off + len - 1;
	    if (len > 40) {        // Big arrays, pseudomedian of 9
		int s = len/8;
		l = med3(x, l,     l+s, l+2*s);
		m = med3(x, m-s,   m,   m+s);
		n = med3(x, n-2*s, n-s, n);
	    }
	    m = med3(x, l, m, n); // Mid-size, med of 3
	}
	double v = x[m];

	// Establish Invariant: v* (<v)* (>v)* v*
	int a = off, b = a, c = off + len - 1, d = c;
	while(true) {
	    while (b <= c && x[b] <= v) {
		if (x[b] == v)
		    swap(x, y, a++, b);
		b++;
	    }
	    while (c >= b && x[c] >= v) {
		if (x[c] == v)
		    swap(x, y, c, d--);
		c--;
	    }
	    if (b > c)
		break;
	    swap(x, y, b++, c--);
	}

	// Swap partition elements back to middle
	int s, n = off + len;
	s = Math.min(a-off, b-a  );  vecswap(x, y, off, b-s, s);
	s = Math.min(d-c,   n-d-1);  vecswap(x, y, b,   n-s, s);

	// Recursively sort non-partition-elements
	if ((s = b-a) > 1)
	    sort1(x, y, off, s);
	if ((s = d-c) > 1)
	    sort1(x, y, n-s, s);
    }

    private static int med3(double x[], int a, int b, int c) {
	return (x[a] < x[b] ?
		(x[b] < x[c] ? b : x[a] < x[c] ? c : a) :
		(x[b] > x[c] ? b : x[a] > x[c] ? c : a));
    }

    private static void swap(double x[], String[] y, int a, int b) {
	double t = x[a];
        String u = y[a];
        y[a] = y[b];
	x[a] = x[b];
        y[b] = u;
	x[b] = t;
    }

    private static void vecswap(double x[], String[] y, int a, int b, int n) {
	for (int i=0; i<n; i++, a++, b++)
	    swap(x, y, a, b);
    }

    private static int binarySearch0(double[] a, int fromIndex, int toIndex,
				     double key) {
	int low = fromIndex;
	int high = toIndex - 1;

	while (low <= high) {
	    int mid = (low + high) >>> 1;
	    double midVal = a[mid];

            int cmp;
            if (midVal < key) {
                cmp = -1;   // Neither val is NaN, thisVal is smaller
            } else if (midVal > key) {
                cmp = 1;    // Neither val is NaN, thisVal is larger
            } else {
                long midBits = Double.doubleToLongBits(midVal);
                long keyBits = Double.doubleToLongBits(key);
                cmp = (midBits == keyBits ?  0 : // Values are equal
                       (midBits < keyBits ? -1 : // (-0.0, 0.0) or (!NaN, NaN)
                        1));                     // (0.0, -0.0) or (NaN, !NaN)
            }

	    if (cmp < 0)
		low = mid + 1;
	    else if (cmp > 0)
		high = mid - 1;
	    else
		return mid; // key found
	}
	return -(low + 1);  // key not found.
    }

}
