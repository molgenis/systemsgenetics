/*
 * GeneLocationObjectSorter.java
 *
 * Created on 23 December 2003, 17:14
 */

package nl.systemsgenetics.eqtlinteractionanalyser.eqtlinteractionanalyser;

/**
 *
 * @author  Like
 */
public class DoubleArrayIntegerObjectSorter extends VectorSorter {
    
    /** Creates a new instance of GeneLocationObjectSorter */
    public DoubleArrayIntegerObjectSorter() {
        super();
    }

    /** Override object comparer
     * @param a the first GeneLocationObject to be compared
     * @param b the second GeneLocationObject to be compared
     * @return true if the first GeneLocationObject.getChrStart() is lower than the second one
     */
    protected boolean lt (Object a, Object b) {
        return (((DoubleArrayIntegerObject)a).intValue < ((DoubleArrayIntegerObject)b).intValue);
    }
    
}