package umcg.genetica.util;

import org.junit.Test;

import java.util.Arrays;
import java.util.Random;

import static org.junit.Assert.assertArrayEquals;

/**
 * Tests for the sorter.
 *
 * @author Kenneth Powers
 */
public class MultiThreadedInPlaceQuickSortTest {
    /**
     * Tests the Sorter implementation on small set of values.
     */
    @Test
    public void testSort() {
        // Define some values to sort.
        Integer[] values = {0, 9, 2, 5, 1, 5, 2, 7, 9, 2, 5, 0, 5, 4, 1};
        // Copy the values into a new array to sort with Java's built in implementation for comparison
        Integer[] sortedValues = new Integer[values.length];
        System.arraycopy(values, 0, sortedValues, 0, values.length);
        // Java sort
        Arrays.sort(sortedValues);
        // Multi-Threaded Quick Sort
        MultiThreadedInPlaceQuickSort.quicksort(values);
        // Assert equality.
        assertArrayEquals(sortedValues, values);
    }

    /**
     * Tests the Sorter implementation on a larger set of values.
     */
    @Test
    public void largeTestSort() {
        // Generate an array of one million random integers.
        Random random = new Random(System.currentTimeMillis());
        Integer[] values = new Integer[10000000];
        final long generateStart = System.currentTimeMillis();
        for (int i = 0; i < values.length; i++) {
            values[i] = random.nextInt();
        }
        final long generateEnd = System.currentTimeMillis();
        System.out.printf("Generation time: %d (%d Integer Objects)%n", generateEnd - generateStart, values.length);

        // Copy the values into a new array to sort with Java's built in implementation for comparison
        Integer[] sortedValues = new Integer[10000000];
        final long copyStart = System.currentTimeMillis();
        System.arraycopy(values, 0, sortedValues, 0, values.length);
        final long copyEnd = System.currentTimeMillis();
        System.out.println("Copy time: " + (copyEnd - copyStart));

        // Start sorting
        final long sortStartTime = System.currentTimeMillis();
        // Java sort
        Arrays.sort(sortedValues);
        final long sortTradeOffTime = System.currentTimeMillis();
        // Multi-Threaded Quick Sort
        MultiThreadedInPlaceQuickSort.quicksort(values);
        final long sortEndTime = System.currentTimeMillis();

        // Assert equality.
        assertArrayEquals(sortedValues, values);

        // Print Statistics
        System.out.println("Java built in implementation: " + (sortTradeOffTime - sortStartTime));
        System.out.println("Multi-threaded quick sort: " + (sortEndTime - sortTradeOffTime));
    }
}
