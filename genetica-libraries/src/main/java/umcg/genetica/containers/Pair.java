/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.containers;

/**
 *
 * @author harm-jan taken from :
 * http://stackoverflow.com/questions/521171/a-java-collection-of-value-pairs-tuples
 */
public class Pair<L, R> implements Comparable<Pair<L, R>> {

    private final L left;
    private final R right;
    private String sep = "-";
    private SORTBY sorter = SORTBY.BOTH;

    public enum SORTBY {

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
        return left.hashCode() ^ right.hashCode();
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
        return this.left.equals(pairo.getLeft())
                && this.right.equals(pairo.getRight());
    }

    @Override
    public int compareTo(Pair<L, R> t) {
        if (this.equals(t)) {
            return 0;
        } else if (sorter == SORTBY.LEFT) {
            if (t.left instanceof Double && this.left instanceof Double) {
                return (Double) t.left > (Double) this.left ? 1 : -1;
            }
            if (t.left instanceof Integer && this.left instanceof Integer) {
                return (Integer) t.left > (Integer) this.left ? 1 : -1;
            }
        } else if (sorter == SORTBY.RIGHT) {
            if (t.right instanceof Double && this.right instanceof Double) {
                return (Double) t.right > (Double) this.right ? 1 : -1;
            }
            if (t.right instanceof Integer && this.right instanceof Integer) {
                return (Integer) t.right > (Integer) this.right ? 1 : -1;
            }
        } else {
            if (t.left instanceof Double && this.left instanceof Double) {
                if ((Double) t.left > (Double) this.left) {
                    return 1;
                } else if ((Double) t.left < (Double) this.left) {
                    return -1;
                } else {
                    return (Double) t.right > (Double) this.right ? 1 : -1;
                }
            }
            
            if (t.left instanceof Integer && this.left instanceof Integer) {
                if ((Integer) t.left > (Integer) this.left) {
                    return 1;
                } else if ((Integer) t.left < (Integer) this.left) {
                    return -1;
                } else {
                    return (Integer) t.right > (Integer) this.right ? 1 : -1;
                }
            }
        }
        
        // for unsupported types, don't sort at all.
        return 0;

    }
}
