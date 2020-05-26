/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.util;
    import java.util.Random;
/**
 *
 * @author MarcJan
 */
public class InplaceArrayQuickSort {

    public static final Random RND = new Random();
 
    private static void swap(Object[] array, int i, int j) {
        Object tmp = array[i];
        array[i] = array[j];
        array[j] = tmp;
    }
 
    private static <E extends Comparable<? super E>> int partition(E[] array, int begin, int end) {
        int index = begin + RND.nextInt(end - begin + 1);
        E pivot = array[index];
        swap(array, index, end);
        for (int i = index = begin; i < end; ++i) {
            if (array[i].compareTo(pivot) <= 0) {
                swap(array, index++, i);
            }
        }
        swap(array, index, end);
        return (index);
    }
 
    private static <E extends Comparable<? super E>> void qsort(E[] array, int begin, int end) {
        if (end > begin) {
            int index = partition(array, begin, end);
            qsort(array, begin, index - 1);
            qsort(array, index + 1, end);
        }
    }
 
    public static <E extends Comparable<? super E>> void sort(E[] array) {
        qsort(array, 0, array.length - 1);
    }
    
    public static <E extends Comparable<? super E>> void sort(E[] array, int fromIndex, int toIndex) {
        qsort(array, fromIndex, toIndex-1);
    }
 
}
