///*
// * To change this template, choose Tools | Templates
// * and open the template in the editor.
// */
//package umcg.genetica.math.matrix2;
//
//import java.util.ArrayList;
//import java.util.HashSet;
//import java.util.LinkedHashMap;
//
///**
// *
// * @author MarcJan
// */
//public class MergeDoubleMatrices {
//
//    /**
//     * Merge a matrix based on shared column identifiers.
//     *
//     * @param matrixI
//     * @param matrixII
//     * @param removeOldMatrix
//     * @return
//     */
//    public static DoubleMatrixDataset<String, String> mergeMatrixBasedOnColumns(DoubleMatrixDataset<String, String> matrixI, DoubleMatrixDataset<String, String> matrixII, boolean removeOldMatrix) {
//        DoubleMatrixDataset<String, String> newMatrix = null;
//
//        matrixI.OrderOnColumns();
//        matrixII.OrderOnColumns();
//
//        if (matrixI.columns() != matrixII.columns()) {
//            HashSet<String> keepColNames1 = new HashSet<String>();
//            keepColNames1.addAll(matrixI.getColObjects());
//            HashSet<String> keepColNames2 = new HashSet<String>();
//            keepColNames2.addAll(matrixII.getColObjects());
//            keepColNames1.retainAll(keepColNames2);
//
//            MatrixHandling.CreatSubsetBasedOnColumns(matrixI, keepColNames1, false);
//            MatrixHandling.CreatSubsetBasedOnColumns(matrixII, keepColNames1, false);
//        }
//
//        if (matrixI.columns() == 0 || matrixII.columns() == 0) {
//            System.out.println("Warning indivduals merging. No shared columns");
//            System.exit(-1);
//        } else if (matrixI.columns() != matrixII.columns()) {
//            System.out.println("Warning indivduals merging. No equal number of columns");
//            System.exit(-1);
//        }
//
//        HashSet<String> keepRowNames1 = new HashSet<String>();
//        keepRowNames1.addAll(matrixI.getRowObjects());
//        keepRowNames1.addAll(matrixII.getRowObjects());
//
//        HashSet<String> keepRowNames2 = new HashSet<String>();
//
//        for (String key : keepRowNames1) {
//            boolean presentMapI = matrixI.hashRows.containsKey(key);
//            boolean presentMapII = matrixII.hashRows.containsKey(key);
//            if (presentMapI ^ presentMapII) {
//                keepRowNames2.add(key);
//            }
//        }
//
//        if (keepRowNames2.size() > 0) {
//            MatrixHandling.CreatSubsetBasedOnRows(matrixI, keepRowNames2, false);
//            MatrixHandling.CreatSubsetBasedOnRows(matrixII, keepRowNames2, false);
//        }
//
//        keepRowNames1 = null;
//        keepRowNames2 = null;
//
//        double[][] newRawData = new double[(matrixII.rows() + matrixI.rows())][matrixI.columns()];
//        LinkedHashMap<String, Integer> newRowMap = new LinkedHashMap<String, Integer>((matrixII.rows() + matrixI.rows()));
//
//        int tmpPos = 0;
//
//        for (int r = 0; r < matrixI.rows(); ++r) {
//            newRowMap.put(matrixI.getRowObjects().get(r), r);
//            for (int s = 0; s < matrixI.columns(); ++s) {
//                newRawData[r][s] = matrixI.getMatrix().get(r, s);
//            }
//            tmpPos++;
//        }
//        for (int r = 0; r < matrixII.rows(); ++r) {
//            newRowMap.put(matrixII.getRowObjects().get(r), r+tmpPos);
//            for (int s = 0; s < matrixII.columns(); ++s) {
//                newRawData[r + tmpPos][s] = matrixII.getMatrix().get(r, s);
//            }
//        }
//
//        newMatrix.setHashCols(matrixI.getHashCols());
//        newMatrix.setHashRows(newRowMap);
//        newMatrix.rawData = newRawData;
//
//        return (newMatrix);
//    }
//
//    /**
//     * Merge a matrix based on row identifiers.
//     *
//     * @param matrixI
//     * @param matrixII
//     * @param removeOldMatrix
//     * @return
//     */
//    public static DoubleMatrixDataset<String, String> mergeMatrixBasedOnRows(DoubleMatrixDataset<String, String> matrixI, DoubleMatrixDataset<String, String> matrixII, boolean removeOldMatrix) {
//        DoubleMatrixDataset<String, String> newMatrix = new DoubleMatrixDataset<String, String>();
//
//        MatrixHandling.OrderOnRows(matrixI);
//        MatrixHandling.OrderOnRows(matrixII);
//
//        if (matrixI.nrRows != matrixII.nrRows) {
//            HashSet<String> keepRowNames1 = new HashSet<String>();
//            keepRowNames1.addAll(matrixI.rowObjects);
//            HashSet<String> keepRowNames2 = new HashSet<String>();
//            keepRowNames2.addAll(matrixII.rowObjects);
//            keepRowNames1.retainAll(keepRowNames2);
//
//            MatrixHandling.FilterCols(matrixI, keepRowNames1);
//            MatrixHandling.FilterCols(matrixII, keepRowNames1);
//        }
//
//        if (matrixI.nrRows == 0 || matrixII.nrRows == 0) {
//            System.out.println("Warning invlaid merging. No shared rows");
//            System.exit(-1);
//        } else if (matrixI.nrRows != matrixII.nrRows) {
//            System.out.println("Warning invlaid merging. No equal number of rows");
//            System.exit(-1);
//        }
//
//        HashSet<String> keepColNames1 = new HashSet<String>();
//        keepColNames1.addAll(matrixI.colObjects);
//        HashSet<String> keepColNames2 = new HashSet<String>();
//        keepColNames2.addAll(matrixII.colObjects);
//        keepColNames1.retainAll(keepColNames2);
//
//        keepColNames2 = new HashSet<String>();
//
//        for (String key : keepColNames1) {
//            boolean presentMapI = matrixI.hashRows.containsKey(key);
//            boolean presentMapII = matrixII.hashRows.containsKey(key);
//            if (!(presentMapI ^ presentMapII)) {
//                keepColNames2.add(key);
//            }
//        }
//        keepColNames1.removeAll(keepColNames2);
//
//        if (keepColNames2.size() > 0) {
//            MatrixHandling.FilterCols(matrixI, keepColNames1);
//            MatrixHandling.FilterCols(matrixII, keepColNames1);
//        }
//
//        keepColNames1 = null;
//        keepColNames2 = null;
//
//
//        double[][] newRawData = new double[(matrixI.nrRows)][(matrixII.nrCols + matrixI.nrCols)];
//        ArrayList<String> newColIds = new ArrayList<String>((matrixII.nrCols + matrixI.nrCols));
//
//        int tmpPos = 0;
//
//        for (int s = 0; s < matrixI.nrCols; ++s) {
//            newColIds.add(matrixI.colObjects.get(s));
//            for (int r = 0; r < matrixI.nrRows; ++r) {
//                newRawData[r][s] = matrixI.rawData[r][s];
//            }
//            tmpPos++;
//        }
//        for (int s = 0; s < matrixII.nrCols; ++s) {
//            newColIds.add(matrixII.colObjects.get(s));
//            for (int r = 0; r < matrixII.nrRows; ++r) {
//                newRawData[r][s + tmpPos] = matrixII.rawData[r][s];
//            }
//        }
//
//        newMatrix.rowObjects = matrixI.rowObjects;
//        newMatrix.nrRows = matrixI.nrRows;
//
//        if (removeOldMatrix) {
//            matrixI = null;
//            matrixII = null;
//        }
//
//        newMatrix.colObjects = newColIds;
//        newMatrix.nrCols = newColIds.size();
//        newMatrix.rawData = newRawData;
//        newMatrix.recalculateHashMaps();
//
//        return (newMatrix);
//    }
//
//    /**
//     * Merge a set of matrices based on row identifiers.
//     *
//     * @param matrixI
//     * @param matrixII
//     * @param removeOldMatrix
//     * @return
//     */
//    public static DoubleMatrixDataset<String, String> combineBasedOnRows(ArrayList<DoubleMatrixDataset> datasets) {
//        DoubleMatrixDataset<String, String> newMatrix = datasets.get(0);
//
//        if (datasets.size() > 1) {
//            for (int i = 1; i < datasets.size(); ++i) {
//                newMatrix = mergeMatrixBasedOnRows(newMatrix, datasets.get(i), false);
//            }
//        }
//
//        return (newMatrix);
//    }
//
//    /**
//     * Merge a set of matrices based on column identifiers.
//     *
//     * @param matrixI
//     * @param matrixII
//     * @param removeOldMatrix
//     * @return
//     */
//    public static DoubleMatrixDataset<String, String> combineBasedOnCols(ArrayList<DoubleMatrixDataset> datasets) {
//        DoubleMatrixDataset<String, String> newMatrix = datasets.get(0);
//
//        if (datasets.size() > 1) {
//            for (int i = 1; i < datasets.size(); ++i) {
//                newMatrix = mergeMatrixBasedOnColumns(newMatrix, datasets.get(i), false);
//            }
//        }
//
//        return (newMatrix);
//    }
//}
