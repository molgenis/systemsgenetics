/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.matrix2;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import umcg.genetica.math.stats.Descriptives;

/**
 *
 * @author MarcJan
 */
public class MatrixTools {
    public static DoubleMatrix2D centerAndScaleColum(DoubleMatrix2D mat) {
        DenseDoubleMatrix2D newData = new DenseDoubleMatrix2D(mat.rows(), mat.columns());
        System.out.println("Standardizing probe mean and standard deviation");
        for (int c = 0; c < mat.columns(); c++) {
            double mean = Descriptives.mean(mat.viewColumn(c).toArray());
            double stdev = Math.sqrt(Descriptives.variance(mat.viewColumn(c).toArray(), mean));
            for (int r = 0; r < mat.rows(); r++) {
                newData.set(r, c, ((mat.getQuick(r, c)-mean)/stdev));
            }
        }

        return(newData);
    }
    public static DoubleMatrix2D centerColum(DoubleMatrix2D mat) {
        DenseDoubleMatrix2D newData = new DenseDoubleMatrix2D(mat.rows(), mat.columns());
        System.out.println("Standardizing probe mean and standard deviation");
        for (int c = 0; c < mat.columns(); c++) {
            double mean = Descriptives.mean(mat.viewColumn(c).toArray());
            for (int r = 0; r < mat.rows(); r++) {
                newData.set(r, c, ((mat.getQuick(r, c)-mean)));
            }
        }

        return(newData);
    }
    
    public static DoubleMatrix2D scaleColum(DoubleMatrix2D mat) {
        DenseDoubleMatrix2D newData = new DenseDoubleMatrix2D(mat.rows(), mat.columns());
        System.out.println("Standardizing probe mean and standard deviation");
        for (int c = 0; c < mat.columns(); c++) {
            double mean = Descriptives.mean(mat.viewColumn(c).toArray());
            double stdev = Math.sqrt(Descriptives.variance(mat.viewColumn(c).toArray(), mean));
            for (int r = 0; r < mat.rows(); r++) {
                newData.set(r, c, ((mat.getQuick(r, c))/stdev));
            }
        }

        return(newData);
    }
}
