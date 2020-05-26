/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.util;

import java.util.Vector;

/**
 *
 * @author harmjan
 */
public class VectorSorter {

    /****************************************************************************
     * A default constructor. It does nothing.
     ****************************************************************************/
    public VectorSorter() {
    }

    /****************************************************************************
     * Sort the given vector.
     * By default it is assumed that the vector contains elements of type String.
     * If not a subclass must be written which overwrites method
     * <tt>lt(Object,Object)</tt>.
     *<P>
     * @param v a vector to be sorted
     ****************************************************************************/
    public void sort(Vector v) {
	quickSort(v, 0, v.size() - 1);
    }

    /****************************************************************************
     * Compare two objects.
     * <P>
     * By default this method works for Strings. It is meant to be overwritten
     * for other objects.
     * <P>
     * @param a the first object to be compared
     * @param b the second object to be compared
     * @return true if the first object is lower than the second one
     ****************************************************************************/
    protected boolean lt(Object a, Object b) {
	return ((String) a).compareTo((String) b) < 0;
    }

    /****************************************************************************
     * The main algorithm.
     ****************************************************************************/
    private void quickSort(Vector v, int lo0, int hi0) {
	int lo = lo0;
	int hi = hi0;
	Object mid;

	if (hi0 > lo0) {
	    // Arbitrarily establishing partition element as the midpoint of
	    // the array.
	    mid = v.elementAt((lo0 + hi0) / 2);

	    // loop through the array until indices cross
	    while (lo <= hi) {
		// find the first element that is greater than or equal to
		// the partition element starting from the left Index.
		while ((lo < hi0) && lt(v.elementAt(lo), mid)) {
		    ++lo;
		}

		// find an element that is smaller than or equal to
		// the partition element starting from the right Index.
		while ((hi > lo0) && lt(mid, v.elementAt(hi))) {
		    --hi;
		}

		// if the indexes have not crossed, swap
		if (lo <= hi) {
		    swap(v, lo, hi);
		    ++lo;
		    --hi;
		}
	    }


	    // If the right index has not reached the left side of array
	    // must now sort the left partition.
	    if (lo0 < hi) {
		quickSort(v, lo0, hi);
	    }

	    // If the left index has not reached the right side of array
	    // must now sort the right partition.
	    if (lo < hi0) {
		quickSort(v, lo, hi0);
	    }
	}
    }

    private static void swap(Vector a, int i, int j) {
	Object T = a.elementAt(i);
	a.setElementAt(a.elementAt(j), i);
	a.setElementAt(T, j);
    }
}
