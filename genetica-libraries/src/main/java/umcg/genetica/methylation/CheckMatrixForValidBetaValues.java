/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.methylation;

import java.util.ArrayList;
import org.apache.commons.collections.primitives.ArrayDoubleList;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author MarcJan
 */
public class CheckMatrixForValidBetaValues {
    
    /**
     * Check if all probe values are actually between 0 and 1. Otherwise set to NA
     *
     * If not kick out sample
     *
     * @param dataset
     * @param replaceCheckedValuesOutOfRange
     * @return
     */
    public static void checkBetaValues(umcg.genetica.math.matrix2.DoubleMatrixDataset<String, String> dataset, boolean replaceCheckedValuesOutOfRange) {

        for (int c = 0; c < dataset.columns(); ++c) {
            for (int r = 0; r < dataset.rows(); ++r) {
                if (!(dataset.getMatrix().get(r, c) >= 0 && dataset.getMatrix().get(r, c) <= 1)) {
                    dataset.getMatrix().set(r, c, Double.NaN);
                }
            }
        }
    }
    
    /**
     * Check if all probe values are actually between 0 and 1. If
     * replaceCheckedValuesOutOfRange is true values that are known to be
     * allowed are replaced. With 0 if between -0.1 and 0 (Due to background
     * correction). With -999 if either 9 or -3.4E38.
     *
     * If not kick out sample
     *
     * @param dataset
     * @param replaceCheckedValuesOutOfRange
     * @return
     */
    public static ArrayList<String> checkMinAndMaxPerSample(DoubleMatrixDataset<String, String> dataset, boolean replaceCheckedValuesOutOfRange) {

        ArrayList<String> columnsToExclude = new ArrayList<String>();

        for (int c = 0; c < dataset.nrCols; ++c) {
            ArrayDoubleList tmp = new ArrayDoubleList();
            for (int r = 0; r < dataset.nrRows; ++r) {
                if (!(dataset.rawData[r][c] >= 0 && dataset.rawData[r][c] <= 1) && dataset.rawData[r][c] != -999) {
                    if (replaceCheckedValuesOutOfRange) {
                        if (dataset.rawData[r][c] >= -0.01 && dataset.rawData[r][c] <= 1) {
                            dataset.rawData[r][c] = 0;
                        } else if (Math.abs(dataset.rawData[r][c] - (-3.4d * Math.pow(10, 38))) < .0000001) {
                            dataset.rawData[r][c] = -999;
                        } else if ((dataset.rawData[r][c] == 9)) {
                            dataset.rawData[r][c] = -999;
                        } else {
                            System.out.println("This shouldn't be reached");
                            System.out.println("This value reached it though: "+ dataset.rawData[r][c]);
                            System.exit(-1);
                        }
                    } else {
                        tmp.add(dataset.rawData[r][c]);
                    }
                }
            }
            if (tmp.size() != 0) {
                if (tmp.size() > 100) {
                    System.out.println("Excluding due to min and max values probe:\t" + c + "\t" + dataset.colObjects.get(c) + "\t" + tmp.size());
                } else {
                    System.out.println("Excluding due to min and max values probe:\t" + c + "\t" + dataset.colObjects.get(c) + "\t" + tmp.size() + "\t" + tmp.toString());
                }
                columnsToExclude.add(dataset.colObjects.get(c));
            }
        }
        return (columnsToExclude);

    }
}
