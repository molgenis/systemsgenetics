/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats;

import static org.testng.Assert.*;

/**
 *
 * @author MarcJan
 */
public class CorrelationNGTest {

    public CorrelationNGTest() {
    }

    @org.testng.annotations.BeforeMethod
    public void setUpMethod() throws Exception {
    }

    /**
     * Test of rankCorrelate method, of class Correlation.
     */
    @org.testng.annotations.Test
    public void testRankCorrelate() {
    }

    /**
     * Test of correlationToZScore method, of class Correlation.
     */
    @org.testng.annotations.Test
    public void testCorrelationToZScore() {
    }

    /**
     * Test of correlate method, of class Correlation.
     */
    @org.testng.annotations.Test
    public void testCorrelate_4args() {
    }

    /**
     * Test of correlate method, of class Correlation.
     */
    @org.testng.annotations.Test
    public void testCorrelate_doubleArr_doubleArr() {
        double[] arrayX = {2.0d, 1.0d, 1.5d, 1.9d, 2.0d};
        double[] arrayY = {4.0d, 3.0d, 3.5d, 3.9d, 4.0d};
        double[] arrayZ = {15.0d, 1.0d, 28.5d, 1.9d, 2.0d};
        double[] arrayA = {81.0d, 83.0d, 0.5d, 53.9d, 6.0d};
        double[] arrayB = {0.0d, 0.0d, 0.0d, 0.0d, 0.0d};
        
        double cor = Correlation.correlate(arrayX, arrayY);

        assertEquals(cor, 1.0d, 0.000001);

        cor = Correlation.correlate(arrayY, arrayX);

        assertEquals(cor, 1.0d, 0.000001);

        cor = Correlation.correlate(arrayY, arrayZ);

        assertEquals(cor, 0.002309785, 0.000001);
        
        cor = Correlation.correlate(arrayZ, arrayA);
        
        assertEquals(cor, -0.3902794, 0.000001);
        
        cor = Correlation.correlate(arrayX, arrayA);
        
        assertEquals(cor, -0.2448009d, 0.000001);
        
        cor = Correlation.correlate(arrayB, arrayA);
        
        assertEquals(Double.isNaN(cor), true);
    }
    
    /**
     * Test of covariate method, of class Correlation.
     */
    @org.testng.annotations.Test
    public void testCovariate_doubleArr_doubleArr() {
        double[] arrayX = {2.0d, 1.0d, 1.5d, 1.9d, 2.0d};
        double[] arrayY = {4.0d, 3.0d, 3.5d, 3.9d, 4.0d};
        double[] arrayZ = {15.0d, 1.0d, 28.5d, 1.9d, 2.0d};
        double[] arrayA = {81.0d, 83.0d, 0.5d, 53.9d, 6.0d};
        double[] arrayB = {0.0d, 0.0d, 0.0d, 0.0d, 0.0d};
        
        double cov = Correlation.covariate(arrayX, arrayY);

        assertEquals(cov, 0.187d, 0.000001);

        cov = Correlation.covariate(arrayY, arrayX);

        assertEquals(cov, 0.187d, 0.000001);

        cov = Correlation.covariate(arrayY, arrayZ);

        assertEquals(cov, 0.012d, 0.000001);
        
        cov = Correlation.covariate(arrayZ, arrayA);
        
        assertEquals(cov, -186.383, 0.000001);
        
        cov = Correlation.covariate(arrayX, arrayA);
        
        assertEquals(cov, -4.208d, 0.000001);
        
        cov = Correlation.covariate(arrayB, arrayA);
        
        assertEquals(cov, 0.0d, 0.000001);
    }

    /**
     * Test of convertCorrelationToZScore method, of class Correlation.
     */
    @org.testng.annotations.Test
    public void testConvertCorrelationToZScore() {
    }

    /**
     * Test of correlate method, of class Correlation.
     */
    @org.testng.annotations.Test
    public void testCorrelate_6args() {
    }
}