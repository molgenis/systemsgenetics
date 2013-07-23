/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats;

import static org.testng.Assert.*;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;
import umcg.genetica.containers.Pair;

/**
 *
 * @author MarcJan
 */
public class DescriptivesNGTest {
    
    public DescriptivesNGTest() {
    }

    @BeforeMethod
    public void setUpMethod() throws Exception {
    }

    /**
     * Test of mean method, of class Descriptives.
     */
    @Test
    public void testMean_floatArr() {
    }

    /**
     * Test of variance method, of class Descriptives.
     */
    @Test
    public void testVariance_floatArr_double() {
    }

    /**
     * Test of lookupSqrt method, of class Descriptives.
     */
    @Test
    public void testLookupSqrt() {
    }

    /**
     * Test of zScoreToPValue method, of class Descriptives.
     */
    @Test
    public void testZScoreToPValue() {
    }

    /**
     * Test of zScore method, of class Descriptives.
     */
    @Test
    public void testZScore() {
    }

    /**
     * Test of mean method, of class Descriptives.
     */
    @Test
    public void testMean_doubleArr() {
        double[] arrayX = {2.0d, 1.0d, 1.5d, 1.9d, 2.0d};
        double[] arrayA = {81.0d, 83.0d, 0.5d, 53.9d, 6.0d};
        double[] arrayB = {0.0d, 0.0d, 0.0d, 0.0d, 0.0d};
        
        double mean = Descriptives.mean(arrayB);

        assertEquals(mean, 0.0d, 0.000001);
        
        mean = Descriptives.mean(arrayX);

        assertEquals(mean, 1.68d, 0.000001);
        
        mean = Descriptives.mean(arrayA);

        assertEquals(mean, 44.88d, 0.000001);
        
    }

    /**
     * Test of mean method, of class Descriptives.
     */
    @Test
    public void testMean_doubleArr_doubleArr() {
        double[] arrayX = {2.0d, 1.0d, 1.5d, 1.9d, 2.0d};
        double[] arrayA = {81.0d, 83.0d, 0.5d, 53.9d, 6.0d};
        double[] arrayB = {0.0d, 0.0d, 0.0d, 0.0d, 0.0d};
        
        Pair<Double, Double> mean = Descriptives.mean(arrayB, arrayA);

        assertEquals(mean.getLeft(), 0.0d, 0.000001);
        assertEquals(mean.getRight(), 44.88d, 0.000001);
        
        mean = Descriptives.mean(arrayB, arrayX);
        
        assertEquals(mean.getLeft(), 0.0d, 0.000001);
        assertEquals(mean.getRight(), 1.68d, 0.000001);
    }

    /**
     * Test of variance method, of class Descriptives.
     */
    @Test
    public void testVariance_doubleArr() {
        double[] arrayX = {2.0d, 1.0d, 1.5d, 1.9d, 2.0d};
        double[] arrayA = {81.0d, 83.0d, 0.5d, 53.9d, 6.0d};
        double[] arrayB = {0.0d, 0.0d, 0.0d, 0.0d, 0.0d};
        
        double mean = Descriptives.mean(arrayB);
        double variance = Descriptives.variance(arrayB, mean);
        
        assertEquals(variance, 0.0d, 0.000001);
        
        mean = Descriptives.mean(arrayX);
        variance = Descriptives.variance(arrayX, mean);

        assertEquals(variance, 0.187d, 0.000001);
        
        mean = Descriptives.mean(arrayA);
        variance = Descriptives.variance(arrayA, mean);

        assertEquals(variance, 1580.097d, 0.000001);
    }

    /**
     * Test of variance method, of class Descriptives.
     */
    @Test
    public void testVariance_doubleArr_double() {
        
        double[] arrayX = {2.0d, 1.0d, 1.5d, 1.9d, 2.0d};
        double[] arrayA = {81.0d, 83.0d, 0.5d, 53.9d, 6.0d};
        double[] arrayB = {0.0d, 0.0d, 0.0d, 0.0d, 0.0d};
        
        double variance = Descriptives.variance(arrayB);

        assertEquals(variance, 0.0d, 0.000001);
        
        variance = Descriptives.variance(arrayX);

        assertEquals(variance, 0.187d, 0.000001);
        
        variance = Descriptives.variance(arrayA);

        assertEquals(variance, 1580.097d, 0.000001);
    }

    /**
     * Test of getZScorePvalueIndex method, of class Descriptives.
     */
    @Test
    public void testGetZScorePvalueIndex() {
    }

    /**
     * Test of convertZscoreToPvalue method, of class Descriptives.
     */
    @Test
    public void testConvertZscoreToPvalue() {
    }

    /**
     * Test of getSqrt method, of class Descriptives.
     */
    @Test
    public void testGetSqrt() {
    }

    /**
     * Test of sum method, of class Descriptives.
     */
    @Test
    public void testSum() {
    }

    /**
     * Test of absSum method, of class Descriptives.
     */
    @Test
    public void testAbsSum() {
    }
}