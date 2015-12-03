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
public class StringIntegerObjectSorter extends VectorSorter {
    
    private java.text.Collator collatorUS = null;

    /** Creates a new instance of GeneLocationObjectSorter */
    public StringIntegerObjectSorter() {
        super();
        collatorUS = java.text.Collator.getInstance(java.util.Locale.US);
    }

    /** Override object comparer
     * @param a the first GeneLocationObject to be compared
     * @param b the second GeneLocationObject to be compared
     * @return true if the first GeneLocationObject.getChrStart() is lower than the second one
     */
    protected boolean lt (Object a, Object b) {
        if (collatorUS.compare(((StringIntegerObject)a).stringValue, ((StringIntegerObject)b).stringValue) >= 0) {
            return false;
        } else {
            return true;
        }
    }
    
}