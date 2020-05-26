/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.util;

import java.util.Arrays;
import java.util.Random;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;
import static org.testng.Assert.*;

/**
 *
 * @author MarcJan
 */
public class MultiThreadedInPlaceQuickSortNGTest {
    
    public MultiThreadedInPlaceQuickSortNGTest() {
    }

    /**
     * Test of quicksort method, of class MultiThreadedInPlaceQuickSort.
     */
    @Test
    public void testQuicksort() {
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
        assertEquals(sortedValues, values);
    }

    public void performanceTestingQuickSorters() {
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
        Integer[] values2 = new Integer[10000000];
        Integer[] values3 = new Integer[10000000];
        final long copyStart = System.currentTimeMillis();
        System.arraycopy(values, 0, values3, 0, values.length);
        System.arraycopy(values, 0, values2, 0, values.length);
        final long copyEnd = System.currentTimeMillis();
        System.out.println("Copy time: " + (copyEnd - copyStart));

        // Start sorting
        final long sortStartTime = System.currentTimeMillis();
        // Java sort
        Arrays.sort(values3);
        final long sortTradeOffTime = System.currentTimeMillis();
        // Multi-Threaded Quick Sort
        MultiThreadedInPlaceQuickSort.quicksort(values);
        final long sortTradeOffTime2 = System.currentTimeMillis();
        // Multi-Threaded Quick Sort
        InplaceArrayQuickSort.sort(values2);
        final long sortEndTime = System.currentTimeMillis();
        // Assert equality.
        assertEquals(values3, values);
        assertEquals(values2, values);

        // Print Statistics
        System.out.println("Java built in implementation: " + (sortTradeOffTime - sortStartTime));
        System.out.println("InplaceArrayQuickSort quick sort: " + (sortEndTime - sortTradeOffTime2));
        System.out.println("Multi-threaded quick sort: " + (sortTradeOffTime2 - sortTradeOffTime));
    }
}