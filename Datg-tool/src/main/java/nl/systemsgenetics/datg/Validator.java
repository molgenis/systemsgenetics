package nl.systemsgenetics.datg;

import org.apache.commons.math3.util.Precision;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.io.IOException;

public class Validator {

    public static void compareTwoMatrices(DoubleMatrixDataset<String, String> m1, DoubleMatrixDataset<String, String> m2) throws IOException {
        compareTwoMatrices(m1, m2, 0.00000001);
    }

    public static void compareTwoMatrices(DoubleMatrixDataset<String, String> m1, DoubleMatrixDataset<String, String> m2, double delta) throws IOException {

        if(m1.rows() != m2.rows()){
            throw new IOException("Rows not equal");
        }
        if(m1.columns()!= m2.columns()){
            throw new IOException("Cols not equal");
        }

        if(!m1.getRowObjects().equals(m2.getRowObjects())){
            throw new IOException("Row names not equal");
        }

        if(!m1.getColObjects().equals(m2.getColObjects())){
            throw new IOException("Col names not equal");
        }

        for (int r = 0; r < m1.rows(); ++r) {
            for (int c = 0; c < m1.columns(); ++c) {
                if(!Precision.equalsIncludingNaN(m1.getElementQuick(r, c), m2.getElementQuick(r, c), delta)){
                    throw new IOException("Difference at r: " + r + " c: " + c);
                }

            }
        }

    }

}
