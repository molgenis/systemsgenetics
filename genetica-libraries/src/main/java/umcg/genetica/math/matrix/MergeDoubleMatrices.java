/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.matrix;

import java.util.ArrayList;
import java.util.HashSet;


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
    public static DoubleMatrixDataset<String, String> mergeMatrixBasedOnColumns(DoubleMatrixDataset<String, String> matrixI, DoubleMatrixDataset<String, String> matrixII, boolean removeOldMatrix) {
        DoubleMatrixDataset<String, String> newMatrix = new DoubleMatrixDataset<String, String>();

        MatrixHandeling.OrderOnColumns(matrixI);
        MatrixHandeling.OrderOnColumns(matrixII);

        if (matrixI.nrCols != matrixII.nrCols) {
            HashSet<String> keepColNames1 = new HashSet<String>();
            keepColNames1.addAll(matrixI.colObjects);
            HashSet<String> keepColNames2 = new HashSet<String>();
            keepColNames2.addAll(matrixII.colObjects);
            keepColNames1.retainAll(keepColNames2);

            MatrixHandeling.FilterCols(matrixI, keepColNames1);
            MatrixHandeling.FilterCols(matrixII, keepColNames1);
        }

        if (matrixI.nrCols == 0 || matrixII.nrCols == 0) {
            System.out.println("Warning indivduals merging. No shared columns");
            System.exit(-1);
        } else if (matrixI.nrCols != matrixII.nrCols) {
            System.out.println("Warning indivduals merging. No equal number of columns");
            System.exit(-1);
        }
        
        HashSet<String> keepRowNames1 = new HashSet<String>();
        keepRowNames1.addAll(matrixI.rowObjects);
        HashSet<String> keepRowNames2 = new HashSet<String>();
        keepRowNames2.addAll(matrixII.rowObjects);
        keepRowNames1.retainAll(keepRowNames2);
        
        keepRowNames2 = new HashSet<String>();
        
        for(String key : keepRowNames1){
            boolean presentMapI = matrixI.hashRows.containsKey(key);
            boolean presentMapII = matrixII.hashRows.containsKey(key);
            if(!(presentMapI ^ presentMapII)){
                keepRowNames2.add(key);
            }
        }
        keepRowNames1.removeAll(keepRowNames2);
        
        if(keepRowNames2.size()>0){
            MatrixHandeling.FilterRows(matrixI, keepRowNames1);
            MatrixHandeling.FilterRows(matrixII, keepRowNames1);
        }
        
        keepRowNames1 = null;
        keepRowNames2 = null;
        
        double[][] newRawData = new double[(matrixII.nrRows + matrixI.nrRows)][matrixI.nrCols];
        ArrayList<String> newRowIds = new ArrayList<String>((matrixII.nrRows + matrixI.nrRows));

        int tmpPos = 0;

        for (int r = 0; r < matrixI.nrRows; ++r) {
            newRowIds.add(matrixI.rowObjects.get(r));
            for (int s = 0; s < matrixI.nrCols; ++s) {
                newRawData[r][s] = matrixI.rawData[r][s];
            }
            tmpPos++;
        }
        for (int r = 0; r < matrixII.nrRows; ++r) {
            newRowIds.add(matrixII.rowObjects.get(r));
            for (int s = 0; s < matrixII.nrCols; ++s) {
                newRawData[r + tmpPos][s] = matrixII.rawData[r][s];
            }
        }

        newMatrix.colObjects = matrixI.colObjects;
        newMatrix.nrCols = matrixI.nrCols;

        if (removeOldMatrix) {
            matrixI = null;
            matrixII = null;
        }

        newMatrix.rowObjects = newRowIds;
        newMatrix.nrRows = newRowIds.size();
        newMatrix.rawData = newRawData;
        newMatrix.recalculateHashMaps();

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
    public static DoubleMatrixDataset<String, String> mergeMatrixBasedOnRows(DoubleMatrixDataset<String, String> matrixI, DoubleMatrixDataset<String, String> matrixII, boolean removeOldMatrix) {
        DoubleMatrixDataset<String, String> newMatrix = new DoubleMatrixDataset<String, String>();

        MatrixHandeling.OrderOnRows(matrixI);
        MatrixHandeling.OrderOnRows(matrixII);

        if (matrixI.nrRows != matrixII.nrRows) {
            HashSet<String> keepRowNames1 = new HashSet<String>();
            keepRowNames1.addAll(matrixI.rowObjects);
            HashSet<String> keepRowNames2 = new HashSet<String>();
            keepRowNames2.addAll(matrixII.rowObjects);
            keepRowNames1.retainAll(keepRowNames2);

            MatrixHandeling.FilterCols(matrixI, keepRowNames1);
            MatrixHandeling.FilterCols(matrixII, keepRowNames1);
        }

        if (matrixI.nrRows == 0 || matrixII.nrRows == 0) {
            System.out.println("Warning invlaid merging. No shared rows");
            System.exit(-1);
        } else if (matrixI.nrRows != matrixII.nrRows) {
            System.out.println("Warning invlaid merging. No equal number of rows");
            System.exit(-1);
        }
        
        HashSet<String> keepColNames1 = new HashSet<String>();
        keepColNames1.addAll(matrixI.colObjects);
        HashSet<String> keepColNames2 = new HashSet<String>();
        keepColNames2.addAll(matrixII.colObjects);
        keepColNames1.retainAll(keepColNames2);
        
        keepColNames2 = new HashSet<String>();
        
        for(String key : keepColNames1){
            boolean presentMapI = matrixI.hashRows.containsKey(key);
            boolean presentMapII = matrixII.hashRows.containsKey(key);
            if(!(presentMapI ^ presentMapII)){
                keepColNames2.add(key);
            }
        }
        keepColNames1.removeAll(keepColNames2);
        
        if(keepColNames2.size()>0){
            MatrixHandeling.FilterCols(matrixI, keepColNames1);
            MatrixHandeling.FilterCols(matrixII, keepColNames1);
        }
        
        keepColNames1 = null;
        keepColNames2 = null;
        
        
        double[][] newRawData = new double[(matrixI.nrRows)][(matrixII.nrCols + matrixI.nrCols)];
        ArrayList<String> newColIds = new ArrayList<String>((matrixII.nrCols + matrixI.nrCols));

        int tmpPos = 0;

        for (int s = 0; s < matrixI.nrCols; ++s) {
            newColIds.add(matrixI.colObjects.get(s));
            for (int r = 0; r < matrixI.nrRows; ++r) {
                newRawData[r][s] = matrixI.rawData[r][s];
            }
            tmpPos++;
        }
        for (int s = 0; s < matrixII.nrCols; ++s) {
            newColIds.add(matrixII.colObjects.get(s));
            for (int r = 0; r < matrixII.nrRows; ++r) {
                newRawData[r][s + tmpPos] = matrixII.rawData[r][s];
            }
        }

        newMatrix.rowObjects = matrixI.rowObjects;
        newMatrix.nrRows = matrixI.nrRows;

        if (removeOldMatrix) {
            matrixI = null;
            matrixII = null;
        }

        newMatrix.colObjects = newColIds;
        newMatrix.nrCols = newColIds.size();
        newMatrix.rawData = newRawData;
        newMatrix.recalculateHashMaps();

        return (newMatrix);
    }

    /**
     * Merge a set of matrices based on row identifiers.
     * 
     * @param matrixI
     * @param matrixII
     * @param removeOldMatrix
     * @return
     */
    public static DoubleMatrixDataset<String, String> combineBasedOnRows(ArrayList<DoubleMatrixDataset> datasets) {
        DoubleMatrixDataset<String, String> newMatrix = datasets.get(0);

        if (datasets.size() > 1) {
            for (int i = 1; i < datasets.size(); ++i) {
                newMatrix = mergeMatrixBasedOnRows(newMatrix, datasets.get(i), false);
            }
        }

        return (newMatrix);
    }

    /**
     * Merge a set of matrices based on column identifiers.
     *
     * @param matrixI
     * @param matrixII
     * @param removeOldMatrix
     * @return
     */
    public static DoubleMatrixDataset<String, String> combineBasedOnCols(ArrayList<DoubleMatrixDataset> datasets) {
        DoubleMatrixDataset<String, String> newMatrix = datasets.get(0);

        if (datasets.size() > 1) {
            for (int i = 1; i < datasets.size(); ++i) {
                newMatrix = mergeMatrixBasedOnColumns(newMatrix, datasets.get(i), false);
            }
        }

        return (newMatrix);
    }
}
