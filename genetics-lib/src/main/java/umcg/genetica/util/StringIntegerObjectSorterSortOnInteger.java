package umcg.genetica.util;

import umcg.genetica.containers.StringIntegerObject;

public class StringIntegerObjectSorterSortOnInteger extends VectorSorter {

    private java.text.Collator collatorUS = null;

    /** Creates a new instance of GeneLocationObjectSorter */
    public StringIntegerObjectSorterSortOnInteger() {
        super();
        collatorUS = java.text.Collator.getInstance(java.util.Locale.US);
    }

    /** Override object comparer
     * @param a the first GeneLocationObject to be compared
     * @param b the second GeneLocationObject to be compared
     * @return true if the first GeneLocationObject.getChrStart() is lower than the second one
     */
    protected boolean lt (Object a, Object b) {
        return ((StringIntegerObject)a).intValue < ((StringIntegerObject)b).intValue;
    }

}