package nl.systemsgenetics.downstreamer.gene;

// TODO: Move this to IO package?
public class IndexedDouble implements Comparable<IndexedDouble> {

    private Double value;
    private int index;

    public IndexedDouble(double value, int index) {
        this.value = value;
        this.index = index;
    }

    public double getValue() {
        return value;
    }

    public void setValue(double value) {
        this.value = value;
    }

    public int getIndex() {
        return index;
    }

    public void setIndex(int index) {
        this.index = index;
    }

    @Override
    public int compareTo(IndexedDouble o) {
        return this.value.compareTo(o.getValue());

    }
}
