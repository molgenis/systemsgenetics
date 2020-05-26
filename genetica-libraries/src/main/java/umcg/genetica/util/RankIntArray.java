/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.util;

/**
 *
 * @author lude
 */
public class RankIntArray {

    public int[] x = null;
    public int[] y = null;
    public cern.colt.Swapper swapper = null;
    public cern.colt.function.tint.IntComparator comp = null;

    public RankIntArray() {
        swapper = new cern.colt.Swapper() {

            public void swap(int a, int b) {
                int t1;
                int t2;
                t1 = x[a];
                x[a] = x[b];
                x[b] = t1;
                t2 = y[a];
                y[a] = y[b];
                y[b] = t2;
            }
        };

        comp = new cern.colt.function.tint.IntComparator() {

            public int compare(int a, int b) {
                return x[a] == x[b] ? 0 : (x[a] < x[b] ? -1 : 1);
            }
        };
    }

    public int[] rank(int[] x) {
        this.x = x.clone();
        y = new int[x.length];
        for (int v = 0; v < x.length; v++) {
            y[v] = v;
        }
        
        cern.colt.GenericSorting.quickSort(0, x.length, comp, swapper);
        
        int[] rank = new int[x.length];
        for (int v = 0; v < x.length; v++) {
            rank[y[v]] = v;
        }
        return rank;
    }
}
