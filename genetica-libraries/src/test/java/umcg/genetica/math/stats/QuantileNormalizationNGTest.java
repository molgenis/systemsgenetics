/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats;

import static org.testng.Assert.*;

import org.testng.annotations.Test;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 * @author MarcJan
 */
public class QuantileNormalizationNGTest {

    public QuantileNormalizationNGTest() {
    }

    /**
     * Test of quantilenormalize method, of class QuantileNormalization.
     */
    @Test
    public void testQuantilenormalize() {

        double[][] testMatrix = {
                {1.0d, 3.0d, 4.0d, 3.0d, 5.0d, 6.0d, 10.0d, 9.0d, 15.0d, 15.0d},
                {0.0d, 4.0d, 1.0d, 15.0d, 8.0d, 6.5d, 6.5d, 6.5d, 1.0d, 0.0d},
                {12.0d, 30.0d, 42.0d, 35.0d, 25.0d, 61.0d, 1.0d, 0.0d, 1.0d, 1.0d},
                {20.0d, 0.0d, 0.0d, 0.0d, 15.0d, 16.0d, 11.0d, 19.0d, 18.0d, 18.0d},
                {12.0d, 30.0d, 42.0d, 35.0d, 25.0d, 61.0d, 1.0d, 0.0d, 1.0d, 1.0d},
                {0.0d, 4.0d, 1.0d, 15.0d, 8.0d, 6.5d, 6.5d, 6.5d, 1.0d, 0.0d}};

        DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<>();
        ds.setMatrix(testMatrix);
        QuantileNormalization.quantilenormalize(ds);

        double[][] expectedMatrix = {
                {5.05d, 2.35d, 8.1d, 2.35d, 1.3d, 1.3d, 25.4d, 25.4d, 25.4d, 25.4d},
                {1.825d, 6.575d, 3.7d, 6.575d, 3.7d, 3.7d, 6.575d, 6.575d, 3.7d, 1.825d},
                {16.75d, 26.65d, 26.65d, 26.65d, 26.65d, 26.65d, 1.825d, 1.825d, 3.7d, 6.575d},
                {27.9d, 1.3d, 1.3d, 1.3d, 8.1d, 8.1d, 27.9d, 27.9d, 27.9d, 27.9d},
                {16.75d, 26.65d, 26.65d, 26.65d, 26.65d, 26.65d, 1.825d, 1.825d, 3.7d, 6.575d},
                {1.825d, 6.575d, 3.7d, 6.575d, 3.7d, 3.7d, 6.575d, 6.575d, 3.7d, 1.825d}};

        for (int i = 0; i < testMatrix.length; ++i) {
            for (int j = 0; j < testMatrix[i].length; ++j) {
                assertEquals(ds.getElement(i, j), expectedMatrix[i][j], 0.0000000001);
            }
            System.out.println("");
        }


    }


}
