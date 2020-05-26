package umcg.genetica.util;

import umcg.genetica.containers.StringDoubleObject;

public class StringDoubleObjectSorterSortOnDouble extends VectorSorter {

    private java.text.Collator collatorUS = null;

    /** Creates a new instance of GeneLocationObjectSorter */
    public StringDoubleObjectSorterSortOnDouble() {
        super();
        collatorUS = java.text.Collator.getInstance(java.util.Locale.US);
    }

    /** Override object comparer
     * @param a the first GeneLocationObject to be compared
     * @param b the second GeneLocationObject to be compared
     * @return true if the first GeneLocationObject.getChrStart() is lower than the second one
     */
    protected boolean lt (Object a, Object b) {
        return ((StringDoubleObject)a).doubleValue < ((StringDoubleObject)b).doubleValue;
    }

}