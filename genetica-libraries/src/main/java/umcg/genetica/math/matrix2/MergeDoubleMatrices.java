/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.matrix2;

import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseLargeDoubleMatrix2D;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author MarcJan
 */
public class MergeDoubleMatrices {

    /**
     * Merge a matrix based on shared column identifiers.
     *
     * @param matrixI
     * @param matrixII
     * @param removeOldMatrix
     * @return
     */
    public static DoubleMatrixDataset<String, String> mergeMatrixBasedOnColumns(DoubleMatrixDataset<String, String> matrixI, DoubleMatrixDataset<String, String> matrixII, boolean removeOldMatrix) throws Exception {
        DoubleMatrixDataset<String, String> newMatrix = null;

        matrixI.OrderOnColumns();
        matrixII.OrderOnColumns();

        if (matrixI.columns() != matrixII.columns()) {
            HashSet<String> keepColNames1 = new HashSet<String>();
            keepColNames1.addAll(matrixI.getColObjects());
            HashSet<String> keepColNames2 = new HashSet<String>();
            keepColNames2.addAll(matrixII.getColObjects());
            keepColNames1.retainAll(keepColNames2);

            matrixI = MatrixHandling.CreatSubsetBasedOnColumns(matrixI, keepColNames1, false);
            matrixII = MatrixHandling.CreatSubsetBasedOnColumns(matrixII, keepColNames1, false);
        }

        if (matrixI.columns() == 0 || matrixII.columns() == 0) {
            System.out.println("Warning indivduals merging. No shared columns");
            throw new Exception("Warning invalid merging. No shared columns");
        } else if (matrixI.columns() != matrixII.columns()) {
            System.out.println("Warning indivduals merging. No equal number of columns");
            throw new Exception("Warning indivduals merging. No equal number of columns");
        }

        HashSet<String> keepRowNames1 = new HashSet<String>();
        keepRowNames1.addAll(matrixI.getRowObjects());
        keepRowNames1.addAll(matrixII.getRowObjects());

        HashSet<String> keepRowNames2 = new HashSet<String>();

        for (String key : keepRowNames1) {
            boolean presentMapI = matrixI.hashRows.containsKey(key);
            boolean presentMapII = matrixII.hashRows.containsKey(key);
            if (presentMapI ^ presentMapII) {
                keepRowNames2.add(key);
            }
        }

        if (keepRowNames2.size() > 0) {
            matrixI = MatrixHandling.CreatSubsetBasedOnRows(matrixI, keepRowNames2, false);
            matrixII = MatrixHandling.CreatSubsetBasedOnRows(matrixII, keepRowNames2, false);
        }

        keepRowNames1 = null;
        keepRowNames2 = null;

        double[][] newRawData = new double[(matrixI.rows()+matrixII.rows())][matrixI.columns()];
        LinkedHashMap<String, Integer> newRowMap = new LinkedHashMap<String, Integer>((matrixI.rows()+matrixII.rows()));

        int tmpPos = 0;

        for (int r = 0; r < matrixI.rows(); ++r) {
            newRowMap.put(matrixI.getRowObjects().get(r), r);
            for (int s = 0; s < matrixI.columns(); ++s) {
                newRawData[r][s] = matrixI.getMatrix().get(r, s);
            }
            tmpPos++;
        }
        for (int r = 0; r < matrixII.rows(); ++r) {
            newRowMap.put(matrixII.getRowObjects().get(r), r+tmpPos);
            for (int s = 0; s < matrixII.columns(); ++s) {
                newRawData[r + tmpPos][s] = matrixII.getMatrix().get(r, s);
            }
        }
        
        
        if ((newRowMap.size() * matrixI.columns()) < (Integer.MAX_VALUE - 2)) {
            newMatrix = new SmallDoubleMatrixDataset<String, String>(new DenseDoubleMatrix2D(newRawData), newRowMap, matrixI.getHashCols());
        } else {
            DenseLargeDoubleMatrix2D matrix = new DenseLargeDoubleMatrix2D(newRowMap.size(), matrixI.columns());
            matrix.assign(newRawData);
            newMatrix = new LargeDoubleMatrixDataset<String, String>(matrix, newRowMap, matrixI.getHashCols());
        }

        return (newMatrix);
    }

    /**
     * Merge a matrix based on row identifiers.
     *
     * @param matrixI
     * @param matrixII
     * @param removeOldMatrix
     * @return
     */
    public static DoubleMatrixDataset<String, String> mergeMatrixBasedOnRows(DoubleMatrixDataset<String, String> matrixI, DoubleMatrixDataset<String, String> matrixII, boolean removeOldMatrix) throws Exception {
        DoubleMatrixDataset<String, String> newMatrix = null;

        matrixI.OrderOnRows();
        matrixII.OrderOnRows();

        if (matrixI.rows() != matrixII.rows()) {
            HashSet<String> keepRowNames1 = new HashSet<String>();
            keepRowNames1.addAll(matrixI.getRowObjects());
            HashSet<String> keepRowNames2 = new HashSet<String>();
            keepRowNames2.addAll(matrixII.getRowObjects());
            keepRowNames1.retainAll(keepRowNames2);
            if(keepRowNames1.size() != matrixI.rows()){
                matrixI = MatrixHandling.CreatSubsetBasedOnRows(matrixI, keepRowNames1, false);
            }
            if(keepRowNames1.size() != matrixII.rows()){
                matrixII = MatrixHandling.CreatSubsetBasedOnRows(matrixII, keepRowNames1, false);
            }
        }

        if (matrixI.rows() == 0 || matrixII.rows() == 0) {
            System.out.println("Warning invalid merging. No shared rows");
            throw new Exception("Warning invalid merging. No shared rows");
        } else if (matrixI.rows() != matrixII.rows()) {
            System.out.println("Warning invalid merging. No equal number of rows");
            throw new Exception("Warning invalid merging. No equal number of rows");
        }

        
        HashSet<String> keepColNames1 = new HashSet<String>();
        keepColNames1.addAll(matrixI.getColObjects());
        keepColNames1.addAll(matrixII.getColObjects());

        HashSet<String> keepColNames2 = new HashSet<String>();

        for (String key : keepColNames1) {
            boolean presentMapI = matrixI.hashRows.containsKey(key);
            boolean presentMapII = matrixII.hashRows.containsKey(key);
            if (presentMapI ^ presentMapII) {
                keepColNames2.add(key);
            }
        }

        if (keepColNames2.size() > 0) {
            matrixI = MatrixHandling.CreatSubsetBasedOnRows(matrixI, keepColNames2, false);
            matrixII = MatrixHandling.CreatSubsetBasedOnRows(matrixII, keepColNames2, false);
        }
        
        keepColNames1 = null;
        keepColNames2 = null;


        double[][] newRawData = new double[(matrixI.rows())][(matrixII.columns() + matrixI.columns())];
        LinkedHashMap<String, Integer> newColMap = new LinkedHashMap<String, Integer>((matrixII.columns() + matrixI.columns()));

        int tmpPos = 0;

        for (int s = 0; s < matrixI.columns(); ++s) {
            newColMap.put(matrixI.getColObjects().get(s), s);
            for (int r = 0; r < matrixI.rows(); ++r) {
                newRawData[r][s] = matrixI.getMatrix().get(r, s);
            }
            tmpPos++;
        }
        for (int s = 0; s < matrixII.columns(); ++s) {
            newColMap.put(matrixII.getColObjects().get(s), s + tmpPos);
            for (int r = 0; r < matrixII.rows(); ++r) {
                newRawData[r][s + tmpPos] = matrixII.getMatrix().get(r, s);
            }
        }

        if ((newColMap.size() * matrixI.columns()) < (Integer.MAX_VALUE - 2)) {
            newMatrix = new SmallDoubleMatrixDataset<String, String>(new DenseDoubleMatrix2D(newRawData), matrixI.getHashRows(), newColMap);
        } else {
            DenseLargeDoubleMatrix2D matrix = new DenseLargeDoubleMatrix2D(matrixI.rows(), newColMap.size());
            matrix.assign(newRawData);
            newMatrix = new LargeDoubleMatrixDataset<String, String>(matrix, matrixI.getHashRows(), newColMap);
        }

        return (newMatrix);
    }

    /**
     * Merge a set of matrices based on row identifiers.
     * Automatic skypping of errors merges.
     *
     * @param matrixI
     * @param matrixII
     * @param removeOldMatrix
     * @return
     */
    public static DoubleMatrixDataset<String, String> combineBasedOnRows(ArrayList<DoubleMatrixDataset<String, String>> datasets) {
        DoubleMatrixDataset<String, String> newMatrix = datasets.get(0);
        System.out.println(newMatrix.getColObjects().size());
        System.out.println(newMatrix.columns());
        
        if (datasets.size() > 1) {
            for (int i = 1; i < datasets.size(); ++i) {
                try {
                    newMatrix = mergeMatrixBasedOnRows(newMatrix, datasets.get(i), false);
                } catch (Exception ex) {
                    Logger.getLogger(MergeDoubleMatrices.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }

        return (newMatrix);
    }

    /**
     * Merge a set of matrices based on column identifiers.
     * Automatic skypping of errors merges.
     *
     * @param matrixI
     * @param matrixII
     * @param removeOldMatrix
     * @return
     */
    public static DoubleMatrixDataset<String, String> combineBasedOnCols(ArrayList<DoubleMatrixDataset<String, String>> datasets) {
        DoubleMatrixDataset<String, String> newMatrix = datasets.get(0);

        if (datasets.size() > 1) {
            for (int i = 1; i < datasets.size(); ++i) {
                try {
                    newMatrix = mergeMatrixBasedOnColumns(newMatrix, datasets.get(i), false);
                } catch (Exception ex) {
                    Logger.getLogger(MergeDoubleMatrices.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }

        return (newMatrix);
    }
}
