/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.matrix2;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseLargeDoubleMatrix2D;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map.Entry;
import org.apache.commons.collections.primitives.ArrayDoubleList;

/**
 *
 * @author MarcJan
 */
public class MatrixHandling {

    /**
     * Remove columns with to many missing values
     *
     * @param dataset DoubleMatrixDataset Expression matrix
     * @param hashColumnsToInclude Ids of samples to include
     * @param missingValue, missing value to check on. If not neccesary put to double.NaN
     */
    public static void RemoveColumnsWithToManyMissingValues(umcg.genetica.math.matrix2.DoubleMatrixDataset<String, String> dataset, int maxMissingValuesPerColumn, double missingValue) {

        HashSet<String> columnsToInclude = new HashSet<String>();

        for (int c = 0; c < dataset.columns(); ++c) {
            int nrMissing = 0;
            for (int r = 0; r < dataset.rows(); ++r) {
                if (dataset.getMatrix().get(r, c) == missingValue || Double.isNaN(dataset.getMatrix().get(r, c))) {
                    nrMissing++;
                }
            }
            if (nrMissing >= maxMissingValuesPerColumn) {
                System.out.println("Excluding:\t" + c + "\t" + dataset.getColObjects().get(c) + "\t" + nrMissing);
            } else {
                columnsToInclude.add(dataset.getColObjects().get(c));
            }
        }

        CreatSubsetBasedOnColumns(dataset, columnsToInclude, true);

    }

    /**
     * Remove rows with to many missing values
     *
     * @param dataset DoubleMatrixDataset Expression matrix
     * @param hashRowsToInclude Ids of rowss to include
     * @param missingValue, missing value to check on. If not neccesary put to double.NaN
     */
    public static void RemoveRowsWithToManyMissingValues(umcg.genetica.math.matrix2.DoubleMatrixDataset<String, String> dataset, int maxMissingValuesPerRow, double missingValue) {

        String[] rowNames = dataset.getRowObjects().toArray(new String[0]);

        HashSet<String> hashRowsToInclude = new HashSet<String>();

        for (int r = 0; r < dataset.rows(); ++r) {
            int nrMissing = 0;
            for (int c = 0; c < dataset.columns(); ++c) {
                if (dataset.getMatrix().get(r, c) == missingValue || Double.isNaN(dataset.getMatrix().get(r, c))) {
                    nrMissing++;
                }
            }
            if (nrMissing >= maxMissingValuesPerRow) {
                System.out.println("Excluding:\t" + r + "\t" + rowNames[r] + "\t" + nrMissing);
            }else {
                hashRowsToInclude.add(rowNames[r]);
            }
        }
        CreatSubsetBasedOnRows(dataset, hashRowsToInclude, false);
    }

    /**
     * Remove identical samples, based on all expression values
     *
     * @param dataset DoubleMatrixDataset Expression matrix
     */
    public static void RemoveDuplicatesSamples(umcg.genetica.math.matrix2.DoubleMatrixDataset<String, String> dataset) {

        HashSet<String> removeEntry = new HashSet<String>();

        for (int c = 0; c < dataset.columns(); ++c) {
            DoubleMatrix1D colc = dataset.getMatrix().viewColumn(c);
            for (int c2 = 0; c2 < dataset.columns(); ++c2) {
                DoubleMatrix1D colc2 = dataset.getMatrix().viewColumn(c2);
                boolean identical = true;
                
                for(int r=0; r < dataset.rows(); ++r){
                    if(colc.getQuick(r) != colc2.getQuick(r)){
                        identical = false;
                        break;
                    }
                }
                
                if(identical){
                    removeEntry.add(dataset.getColObjects().get(c));
                }
            }
        }

        if (removeEntry.size() > 0) {
            RemoveColumns(dataset, removeEntry);
        }
    }
    
    /**
     * Append a static prefix to the column names
     *
     * @param in
     * @param prefix
     */
    public static void appendPrefixToColnames(umcg.genetica.math.matrix2.DoubleMatrixDataset<String, String> in, String prefix) {
        LinkedHashMap<String, Integer> newColObjects = new LinkedHashMap<String, Integer>();
        for (Entry<String, Integer> t : in.getHashCols().entrySet()) {
            StringBuilder colName = new StringBuilder();
            colName.append(prefix);
            colName.append("_");
            colName.append(t);
            newColObjects.put(colName.toString(), t.getValue());
        }
        in.setHashCols(newColObjects);
    }
    
    /**
     * Replace missing values in the double matrix per sample. Using either the
     * mean if useMedian is false or the median is useMedian is true.
     *
     * @param rawData
     * @param useMedian
     * @param NaValue
     */
    public static void ReplaceMissingValuesPerColumn(DoubleMatrix2D rawData, boolean useMedian, double NaValue) {

        for (int s = 0; s < rawData.columns(); ++s) {

            System.out.println("Processing sample: " + s);
            boolean needsReplacement = false;
            ArrayDoubleList nonNAvalues = new ArrayDoubleList();
            
            for (int p = 0; p < rawData.rows(); ++p) {
                if (rawData.get(p, s) == NaValue) {
                    needsReplacement = true;
                } else {
                    nonNAvalues.add(rawData.get(p, s));
                }
            }

            if (needsReplacement) {
                double replacementValue;
                if (useMedian) {
                    replacementValue = JSci.maths.ArrayMath.median(nonNAvalues.toArray(new double[0]));
                } else {
                    replacementValue = JSci.maths.ArrayMath.mean(nonNAvalues.toArray(new double[0]));
                }

                for (int p = 0; p < rawData.rows(); ++p) {
                    if (rawData.get(p, s) == NaValue) {
                        rawData.set(p, s, replacementValue);
                    }
                }
            }
        }
    }
    
    /**
     * Replace missing values in the double matrix per sample. Using either the
     * mean if useMedian is false or the median is useMedian is true.
     *
     * @param rawData
     * @param useMedian
     * @param NaValue
     */
    public static void ReplaceMissingValuesPerRow(DoubleMatrix2D rawData, boolean useMedian, double NaValue) {

        for (int p = 0; p < rawData.rows(); ++p) {

            System.out.println("Processing row: " + p);
            boolean needsReplacement = false;
            ArrayDoubleList nonNAvalues = new ArrayDoubleList();
            
            for (int s = 0; s < rawData.rows(); ++s) {
                if (rawData.get(p, s) == NaValue) {
                    needsReplacement = true;
                } else {
                    nonNAvalues.add(rawData.get(p, s));
                }
            }

            if (needsReplacement) {
                double replacementValue;
                if (useMedian) {
                    replacementValue = JSci.maths.ArrayMath.median(nonNAvalues.toArray(new double[0]));
                } else {
                    replacementValue = JSci.maths.ArrayMath.mean(nonNAvalues.toArray(new double[0]));
                }

                for (int s = 0; s < rawData.rows(); ++s) {
                    if (rawData.get(p, s) == NaValue) {
                        rawData.set(p, s, replacementValue);
                    }
                }
            }
        }
    }
    
    //Hier moeten we nog even aan sleutelen.
    
    /**
     * Remove probes without correct mapping known on forehand
     *
     * @param dataset DoubleMatrixDataset containing the matrix of interest
     * @param probesToBeRemoved ArrayList<String> with identifiers of probes
     * that should be removed
     */
    public static DoubleMatrixDataset<String, String> RemoveProbes(umcg.genetica.math.matrix2.DoubleMatrixDataset<String, String> dataset, HashSet<String> probesToBeRemoved) {
        return CreatSubsetBasedOnRows(dataset, probesToBeRemoved, true);
    }

    public static void fixLinkedHashes(LinkedHashMap<String, Integer> hashMap) {
        int i=0;
        for(Entry<String, Integer> e : hashMap.entrySet()){
           e.setValue(i);
           i++;
        }
    }

    public static void RenameRows(DoubleMatrixDataset<String, ?> dataset, HashMap<String, String> mappedProbeList) {
        LinkedHashMap<String, Integer> newRowNames = new LinkedHashMap<String, Integer>(dataset.rows());
        
        for (Entry<String, Integer> e : dataset.getHashRows().entrySet()) {
            if ((mappedProbeList.containsKey(e.getKey()))) {
                newRowNames.put(mappedProbeList.get(e.getKey()), e.getValue());
            } else {
                newRowNames.put(e.getKey(), e.getValue());
            }
        }

        dataset.setHashRows(newRowNames);
    }

    public static DoubleMatrixDataset<String, String> MergeMatrixWithNonOverlappingColNames(DoubleMatrixDataset<String, String> m1, DoubleMatrixDataset<String, String> m2) {
        HashSet<String> colNames = new HashSet<String>();
        HashSet<String> rowNames = new HashSet<String>();
        
        colNames.addAll(m1.getColObjects());
        colNames.addAll(m2.getColObjects());
        
        rowNames.addAll(m1.getRowObjects());
        rowNames.addAll(m2.getRowObjects());
        
        DoubleMatrixDataset<String, String> mergedMatrix = new DoubleMatrixDataset<String, String>(rowNames, colNames);
        
        for(String colsM1 : m1.getColObjects() ){
            int originalColId = m1.getHashCols().get(colsM1);
            int newColId = mergedMatrix.getHashCols().get(colsM1);
            for(String rowsM1 : m1.getRowObjects()){
                int originalRowId = m1.getHashRows().get(rowsM1);
                int newRowId = mergedMatrix.getHashRows().get(rowsM1);
                mergedMatrix.matrix.setQuick(newRowId, newColId, m1.matrix.getQuick(originalRowId, originalColId));
            }
        }
        
        for(String colsM2 : m2.getColObjects() ){
            int originalColId = m2.getHashCols().get(colsM2);
            int newColId = mergedMatrix.getHashCols().get(colsM2);
            for(String rowsM2 : m2.getRowObjects()){
                int originalRowId = m2.getHashRows().get(rowsM2);
                int newRowId = mergedMatrix.getHashRows().get(rowsM2);
                mergedMatrix.matrix.setQuick(newRowId, newColId, m2.matrix.getQuick(originalRowId, originalColId));
            }
        }
        
        return mergedMatrix;
    }
    
    /**
     * Remove rows without correct mapping known on forehand
     *
     * @param dataset DoubleMatrixDataset containing the matrix of interest
     * @param probesToBeRemoved ArrayList<String> with identifiers of probes
     * that should be removed
     */
    public DoubleMatrixDataset<String, String> RemoveRows(umcg.genetica.math.matrix2.DoubleMatrixDataset<String, String> dataset, HashSet<String> probesToBeRemoved) {
        return CreatSubsetBasedOnRows(dataset, probesToBeRemoved, true);
    }
    
    /**
     * Remove or filter rows.
     *
     * @param dataset DoubleMatrixDataset containing the matrix of interest
     * @param rowNames ArrayList<String> with identifiers of probes that
     * @param removeRows, if true row ids in rowNames are removed if false rowNames are selected others are removed.
     * should be removed
     */
    public static DoubleMatrixDataset<String, String> CreatSubsetBasedOnRows(umcg.genetica.math.matrix2.DoubleMatrixDataset<String, String> dataset, HashSet<String> rowNames, boolean removeRows) {
        LinkedHashMap<String, Integer> rowMap = new LinkedHashMap<String, Integer>();
        
        int newCounter = 0;
        for (String t : dataset.getRowObjects()) {
            if(removeRows && !rowNames.contains(t)){
                rowMap.put(t, newCounter);
                newCounter++;
            } else if(!removeRows && rowNames.contains(t)){
                rowMap.put(t, newCounter);
                newCounter++;
            }
        }
        
        DoubleMatrixDataset<String, String> matrix = new DoubleMatrixDataset<String, String>(rowMap, dataset.getHashCols());
        for (int p = 0; p < matrix.rows(); ++p) {
            int originalRowNumer = dataset.getHashRows().get(matrix.getRowObjects().get(p));
            for (int s = 0; s < matrix.columns(); ++s) {
                matrix.getMatrix().setQuick(p, s, dataset.getMatrix().getQuick(originalRowNumer, s));
            }
        }

        return matrix;
    }

    /**
     * Remove samples on forehand
     *
     * @param dataset DoubleMatrixDataset containing the matrix of interest
     * @param samplesToBeRemoved ArrayList<String> with identifiers of probes
     * that should be removed
     */
    public static DoubleMatrixDataset<String, String> RemoveSamples(umcg.genetica.math.matrix2.DoubleMatrixDataset<String, String> dataset, HashSet<String> samplesToBeRemoved) {
        return CreatSubsetBasedOnColumns(dataset, samplesToBeRemoved, true);
    }
    
    /**
     * Remove samples on forehand
     *
     * @param dataset DoubleMatrixDataset containing the matrix of interest
     * @param samplesToBeRemoved ArrayList<String> with identifiers of probes
     * that should be removed
     */
    public static DoubleMatrixDataset<String, String> RemoveColumns(umcg.genetica.math.matrix2.DoubleMatrixDataset<String, String> dataset, HashSet<String> samplesToBeRemoved) {
        return CreatSubsetBasedOnColumns(dataset, samplesToBeRemoved, true);
    }
    
    /**
     * Filter out columns.
     * Keep all columns that in in the hashset.
     *
     * @param dataset
     * @param colNames
     * @param remove, if true col ids in colNames are removed if false colNames are selected others are removed.
     */
    public static DoubleMatrixDataset<String, String> CreatSubsetBasedOnColumns(umcg.genetica.math.matrix2.DoubleMatrixDataset<String, String> dataset, HashSet<String> colNames, boolean remove) {
        LinkedHashMap<String, Integer> columnMap = new LinkedHashMap<String, Integer>();
        
        int newCounter = 0;
        for (String t : dataset.getColObjects()) {
            if(remove && !colNames.contains(t)){
                columnMap.put(t, newCounter);
                newCounter++;
            } else if(!remove && colNames.contains(t)){
                columnMap.put(t, newCounter);
                newCounter++;
            }
        }
        
        DoubleMatrixDataset<String, String> matrix = new DoubleMatrixDataset<String, String>(dataset.getHashRows(), columnMap);
        for (int s = 0; s < matrix.columns(); ++s) {
            int originalSampleNumer = dataset.getHashCols().get(matrix.getColObjects().get(s));
            for (int p = 0; p < matrix.rows(); ++p) {
                matrix.getMatrix().setQuick(p, s, dataset.getMatrix().getQuick(p, originalSampleNumer));
            }
        }

        return matrix;
    }
    
    public static void RenameCols(umcg.genetica.math.matrix2.DoubleMatrixDataset<?, String> dataset, HashMap<String, String> newNames) {
        LinkedHashMap<String, Integer> newColNames = new LinkedHashMap<String, Integer>(dataset.columns());
        
        for (Entry<String, Integer> e : dataset.getHashCols().entrySet()) {
            if ((newNames.containsKey(e.getKey()))) {
                newColNames.put(newNames.get(e.getKey()), e.getValue());
            } else {
                newColNames.put(e.getKey(), e.getValue());
            }
        }

        dataset.setHashCols(newColNames);
    }
    
    /**
     * Replace zero's with nulls.
     *
     * @param rawData
     */
    public static void ReplaceZerosToNull(DoubleMatrix2D rawData) {
        for (int c = 0; c < rawData.columns(); ++c) {
            for (int r = 0; r < rawData.rows(); ++r) {
                if (rawData.getQuick(r, c) == 0.0d) {
                    rawData.setQuick(r, c, Double.NaN);
                }
            }
        }
    }
    
}
