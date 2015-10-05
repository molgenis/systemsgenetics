/*
 * VectorSorter.java
 *
 * Created on 23 December 2003, 17:02
 */

package nl.systemsgenetics.eqtlinteractionanalyser.eqtlinteractionanalyser;

    import java.util.*;
    import java.io.*;


/** An implementation of a quick sort algorithm.
 * 
 *  Shamelessly taken from Martin Senger's (senger@ebi.ac.uk) Java collection
 *
 *  It sorts a vector of strings, but can be easily extended to sort
 *  a vector of other objects (by overwritting just one method
 *  <b>lt(Object,Object)</b>).
 * <P>
 *  It implements a generic version of C.A.R Hoare's Quick Sort
 *  algorithm.
 * <P>
 *  The code is based on example given in java swing package.
 * <P>
 *
 * This is an example how to use this class to sort a vector of strings:
 * <PRE>
 *    Vector v = new Vector();
 *    v.addElement ("X");
 *    v.addElement ("A");
 *    new Sorter().sort (v);
 * </PRE>
 *
 * @author <A HREF="mailto:senger@ebi.ac.uk">Martin Senger</A>
 * @version $Id: VectorSorter.java,v 1.1.1.1 2004/01/26 09:27:02 lude Exp $
 * @see TestSorter
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
    public void sort (Vector v) {
	quickSort (v, 0, v.size() - 1);
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
    protected boolean lt (Object a, Object b) {
	return ((String)a).compareTo ((String)b) < 0;
    }

    /****************************************************************************
     * The main algorithm.
     ****************************************************************************/
    private void quickSort (Vector v, int lo0, int hi0) {
	int lo = lo0;
	int hi = hi0;
	Object mid;

	if (hi0 > lo0) {
	    // Arbitrarily establishing partition element as the midpoint of
	    // the array.
	    mid = v.elementAt ((lo0 + hi0) / 2);

	    // loop through the array until indices cross
	    while (lo <= hi) {
		// find the first element that is greater than or equal to
		// the partition element starting from the left Index.
		while ((lo < hi0) && lt (v.elementAt (lo), mid)) {
		    ++lo;
		}

		// find an element that is smaller than or equal to
		// the partition element starting from the right Index.
		while ((hi > lo0) && lt (mid, v.elementAt(hi))) {
		    --hi;
		}

		// if the indexes have not crossed, swap
		if (lo <= hi) {
		    swap (v, lo, hi);
		    ++lo;
		    --hi;
		}
	    }


	    // If the right index has not reached the left side of array
	    // must now sort the left partition.
	    if (lo0 < hi) {
		quickSort (v, lo0, hi);
	    }

	    // If the left index has not reached the right side of array
	    // must now sort the right partition.
	    if (lo < hi0) {
		quickSort (v, lo, hi0);
	    }
	}
    }

    private static void swap (Vector a, int i, int j) {
	Object T = a.elementAt(i);
	a.setElementAt (a.elementAt(j), i);
	a.setElementAt (T, j);
    }

}
