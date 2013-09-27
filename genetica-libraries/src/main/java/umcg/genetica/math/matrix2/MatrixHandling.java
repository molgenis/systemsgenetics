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

    private static void fixOrdering(LinkedHashMap<String, Integer> hashMap) {
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
        
        int newSize = 0;
        HashSet<String> removeList = new HashSet<String>();
        
        if(removeRows){
            for (String t : dataset.getRowObjects()) {
                if (!rowNames.contains(t)) {
                    newSize++;
                } else {
                    removeList.add(t);
                }
            }
            if(removeList.isEmpty()){
                return(dataset);
            }
        } else {
            for (String t : dataset.getRowObjects()) {
                if (rowNames.contains(t)) {
                    newSize++;
                } else {
                    removeList.add(t);
                }
            }
            if(newSize == dataset.rows()){
                return(dataset);
            }
        }
        
        DoubleMatrix2D matrix;
        if ((dataset.columns() * (long)newSize) < (Integer.MAX_VALUE - 2)) {
            matrix = new DenseDoubleMatrix2D(newSize, dataset.columns());
        } else{
            matrix = new DenseLargeDoubleMatrix2D(newSize, dataset.columns());
        }
        
        int probeId = -1;
        ArrayList<String> rowObj = dataset.getRowObjects();

        if(removeRows){
            for (int p = 0; p < dataset.rows(); ++p) {
                if (!(rowNames.contains(rowObj.get(p)))) {
                    probeId++;
                    for (int s = 0; s < dataset.columns(); ++s) {
                        matrix.setQuick(p, s, dataset.getMatrix().getQuick(p, s));
                    }
                }
            }
        } else {
            for (int p = 0; p < dataset.rows(); ++p) {
                if ((rowNames.contains(rowObj.get(p)))) {
                    probeId++;
                    for (int s = 0; s < dataset.columns(); ++s) {
                        matrix.setQuick(p, s, dataset.getMatrix().getQuick(p, s));
                    }
                }
            }
        }

        for (String r : removeList) {
            dataset.hashRows.remove(r);
        }

        fixOrdering(dataset.hashRows);
        
        return new DoubleMatrixDataset<String, String>(matrix, dataset.hashRows, dataset.hashCols);
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
        int newSize = 0;
        HashSet<String> removeList = new HashSet<String>();
        
        if(remove){
            for (String t : dataset.getColObjects()) {
                if (!colNames.contains(t)) {
                    newSize++;
                } else {
                    removeList.add(t);
                }   
            }
            if(removeList.isEmpty()){
                return(dataset);
            }
        } else {
            for (String t : dataset.getColObjects()) {
                if (colNames.contains(t)) {
                    newSize++;
                } else {
                    removeList.add(t);
                }
            }
            if(newSize == dataset.columns()){
                return(dataset);
            }
        }
        
        DoubleMatrix2D matrix;
        if ((dataset.rows() * (long)newSize) < (Integer.MAX_VALUE - 2)) {
            matrix = new DenseDoubleMatrix2D(dataset.rows(), newSize);
        } else {
             matrix = new DenseLargeDoubleMatrix2D(dataset.rows(), newSize);
        }
        
        int sampleId = -1;

        ArrayList<String> colObj = dataset.getColObjects();
        if(remove){
            for (int s = 0; s < dataset.columns(); ++s) {
                if (!(colNames.contains(colObj.get(s)))) {
                    sampleId++;
                    for (int p = 0; p < dataset.rows(); ++p) {
                        matrix.setQuick(p, sampleId, dataset.getMatrix().getQuick(p, s));
                    }
                }
            }
        } else {
            for (int s = 0; s < dataset.columns(); ++s) {
                if ((colNames.contains(colObj.get(s)))) {
                    sampleId++;
                    for (int p = 0; p < dataset.rows(); ++p) {
                        matrix.setQuick(p, sampleId, dataset.getMatrix().getQuick(p, s));
                    }
                }
            }
        }
        for (String r : removeList) {
            dataset.hashCols.remove(r);
        }

        fixOrdering(dataset.hashCols);
        
        return new DoubleMatrixDataset<String, String>(matrix, dataset.hashRows, dataset.hashCols);
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
}
