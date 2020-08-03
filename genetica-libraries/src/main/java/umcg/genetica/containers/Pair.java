/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.containers;

/**
 * @author harm-jan taken from :
 * http://stackoverflow.com/questions/521171/a-java-collection-of-value-pairs-tuples
 */
public class Pair<L, R> implements Comparable<Pair<L, R>> {

    private final L left;
    private final R right;
    private String sep = "-";
    private SORTBY sorter = SORTBY.BOTH;

    public static enum SORTBY {

        LEFT, RIGHT, BOTH;
    }

    public Pair(L left, R right) {
        this.left = left;
        this.right = right;
    }

    public Pair(L left, R right, SORTBY sorter) {
        this.left = left;
        this.right = right;
        this.sorter = sorter;
    }

    public Pair(L left, R right, String sep) {
        this(left, right);
        this.sep = sep;
    }

    public Pair(L left, R right, String sep, SORTBY sorter) {
        this.left = left;
        this.right = right;
        this.sorter = sorter;
        this.sep = sep;
    }

    public void setSorter(SORTBY s) {
        this.sorter = s;
    }

    public L getLeft() {
        return left;
    }

    public R getRight() {
        return right;
    }

    @Override
    public String toString() {
        return left.toString() + sep + right.toString();
    }

    @Override
    public int hashCode() {
        if (sorter.equals(SORTBY.BOTH)) {
            return left.hashCode() ^ right.hashCode();
        } else if (sorter.equals(SORTBY.LEFT)) {
            return this.left.hashCode();
        } else {
            return this.right.hashCode();
        }
    }

    @Override
    public boolean equals(Object o) {
        if (o == null) {
            return false;
        }
        if (!(o instanceof Pair)) {
            return false;
        }
        Pair pairo = (Pair) o;
        if (sorter.equals(SORTBY.BOTH)) {
            return this.left.equals(pairo.getLeft())
                    && this.right.equals(pairo.getRight());
        } else if (sorter.equals(SORTBY.LEFT)) {
            return this.left.equals(pairo.getLeft());
        } else {
            return this.right.equals(pairo.getRight());
        }

    }

    @Override
    public int compareTo(Pair<L, R> toCompare) {
        if (this.equals(toCompare)) {
            return 0;
        } else if (sorter.equals(SORTBY.LEFT)) {
            if (toCompare.left instanceof Double && this.left instanceof Double) {
                return ((Double) this.left).compareTo((Double) toCompare.left);
            }
            if (toCompare.left instanceof Integer && this.left instanceof Integer) {
                return ((Integer) this.left).compareTo((Integer) toCompare.left);
            }
        } else if (sorter.equals(SORTBY.RIGHT)) {
            if (toCompare.right instanceof Double && this.right instanceof Double) {
                return ((Double) this.right).compareTo((Double) toCompare.right);
            }
            if (toCompare.right instanceof Integer && this.right instanceof Integer) {
                return ((Integer) this.right).compareTo((Integer) toCompare.right);
            }
        } else {
            if (toCompare.left instanceof Double && this.left instanceof Double) {
                int comp = ((Double) this.left).compareTo((Double) toCompare.left);
                if (comp == 0) {
                    if (toCompare.right instanceof Double && this.right instanceof Double) {
                        return ((Double) this.right).compareTo((Double) toCompare.right);
                    }
                    if (toCompare.right instanceof Integer && this.right instanceof Integer) {
                        return ((Integer) this.right).compareTo((Integer) toCompare.right);
                    }
                } else {
                    return comp;
                }
            }
            if (toCompare.left instanceof Integer && this.left instanceof Integer) {
                int comp = ((Integer) this.left).compareTo((Integer) toCompare.left);
                if (comp == 0) {
                    if (toCompare.right instanceof Double && this.right instanceof Double) {
                        return ((Double) this.right).compareTo((Double) toCompare.right);
                    }
                    if (toCompare.right instanceof Integer && this.right instanceof Integer) {
                        return ((Integer) this.right).compareTo((Integer) toCompare.right);
                    }
                } else {
                    return comp;
                }
            }
        }

        // for unsupported types, don't sort at all.
        return 0;
    }
}
