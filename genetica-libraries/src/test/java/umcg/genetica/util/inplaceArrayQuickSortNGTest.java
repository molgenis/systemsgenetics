/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.util;

import java.util.Arrays;
import static org.testng.Assert.*;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

/**
 *
 * @author MarcJan
 */
public class inplaceArrayQuickSortNGTest {
    
    public inplaceArrayQuickSortNGTest() {
    }

    @BeforeMethod
    public void setUpMethod() throws Exception {
    }

    /**
     * Test of sort method, of class InplaceArrayQuickSort.
     */
    @Test
    public void testSort() {

        Integer[] l1 = { 5, 1024, 1, 88, 0, 1024 };
        Integer[] l1c = { 5, 1024, 1, 88, 0, 1024 };
        
        InplaceArrayQuickSort.sort(l1);
        Arrays.sort(l1c);
        
        boolean correct=true;
        for(int i=0;i<l1.length;i++){
            if(!l1[i].equals(l1c[i])){
                correct = false;
                break;
            }
        }
        assertEquals(correct, true);
        
        Double[] l3 = { 5.0d, 1024.5d, 1.0d, 88.0d, 0.0d, 1024.0d };
        Double[] l3c = { 5.0d, 1024.5d, 1.0d, 88.0d, 0.0d, 1024.0d };
        
        InplaceArrayQuickSort.sort(l3,1,5);
        Arrays.sort(l3c,1,5);
        
        correct=true;
        for(int i=0;i<l3.length;i++){
            if(!l3[i].equals(l3c[i])){
                correct = false;
                break;
            }
        }
        assertEquals(correct, true);
        
        
        
 
        String[] l2 = { "gamma", "beta", "alpha", "zoolander" };
        String[] l2c = { "gamma", "beta", "alpha", "zoolander" };

        InplaceArrayQuickSort.sort(l2);
        Arrays.sort(l2c);
        
        correct=true;
        for(int i=0;i<l2.length;i++){
            if(!l2[i].equals(l2c[i])){
                correct = false;
                break;
            }
        }
        assertEquals(correct, true);
        
    }
}