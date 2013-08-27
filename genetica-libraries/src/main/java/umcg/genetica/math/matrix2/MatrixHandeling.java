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
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map.Entry;
import org.apache.commons.collections.primitives.ArrayDoubleList;

/**
 *
 * @author MarcJan
 */
public class MatrixHandeling {

    /**
     * Remove columns with to many missing values
     *
     * @param dataset DoubleMatrixDataset Expression matrix
     * @param hashColumnsToInclude Ids of samples to include
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

        FilterCols(dataset, columnsToInclude);

    }

    /**
     * Remove rows with to many missing values
     *
     * @param dataset DoubleMatrixDataset Expression matrix
     * @param hashRowsToInclude Ids of rowss to include
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
        FilterRows(dataset, hashRowsToInclude);
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
     * Remove rows without correct mapping known on forehand
     *
     * @param dataset DoubleMatrixDataset containing the matrix of interest
     * @param probesToBeRemoved ArrayList<String> with identifiers of probes
     * that should be removed
     */
    public static DoubleMatrixDataset<String, String> RemoveRows(umcg.genetica.math.matrix2.DoubleMatrixDataset<String, String> dataset, HashSet<String> probesToBeRemoved) {
        int newSize = 0;
        HashSet<String> removeList = new HashSet<String>();

        for (String t : dataset.getRowObjects()) {
            if (!probesToBeRemoved.contains(t)) {
                newSize++;
            } else {
                removeList.add(t);
            }
        }

        double[][] newRawData = new double[newSize][dataset.columns()];

        int probeId = -1;
        ArrayList<String> rowObj = dataset.getRowObjects();

        for (int p = 0; p < dataset.rows(); ++p) {
            if (!(probesToBeRemoved.contains(rowObj.get(p)))) {
                probeId++;
                for (int s = 0; s < dataset.columns(); ++s) {
                    newRawData[probeId][s] = dataset.getMatrix().get(p, s);
                }
            }
        }

        for (String r : removeList) {
            dataset.hashRows.remove(r);
        }

        if ((dataset.columns() * newRawData.length) < (Integer.MAX_VALUE - 2)) {
            dataset = new SmallDoubleMatrixDataset<String, String>(new DenseDoubleMatrix2D(newRawData), dataset.hashRows, dataset.hashCols);
        } else {
            DenseLargeDoubleMatrix2D matrix = new DenseLargeDoubleMatrix2D(newRawData.length, dataset.columns());
            matrix.assign(newRawData);
            dataset = new LargeDoubleMatrixDataset<String, String>(matrix, dataset.hashRows, dataset.hashCols);
        }
        return(dataset);
    }

    /**
     * Remove probes without correct mapping known on forehand
     *
     * @param dataset DoubleMatrixDataset containing the matrix of interest
     * @param probesToBeRemoved ArrayList<String> with identifiers of probes
     * that should be removed
     */
    public static void RemoveProbes(umcg.genetica.math.matrix2.DoubleMatrixDataset<String, String> dataset, HashSet<String> probesToBeRemoved) {
        RemoveRows(dataset, probesToBeRemoved);
    }

    
    /**
     * Remove samples on forehand
     *
     * @param dataset DoubleMatrixDataset containing the matrix of interest
     * @param samplesToBeRemoved ArrayList<String> with identifiers of probes
     * that should be removed
     */
    public static DoubleMatrixDataset<String, String> RemoveColumns(umcg.genetica.math.matrix2.DoubleMatrixDataset<String, String> dataset, HashSet<String> samplesToBeRemoved) {
        
        int newSize = 0;
        HashSet<String> removeList = new HashSet<String>();

        for (String t : dataset.getColObjects()) {
            if (!samplesToBeRemoved.contains(t)) {
                newSize++;
            } else {
                removeList.add(t);
            }
        }

        double[][] newRawData = new double[dataset.rows()][newSize];

        int sampleId = -1;

        ArrayList<String> colObj = dataset.getColObjects();

        for (int s = 0; s < dataset.columns(); ++s) {
            if (!(samplesToBeRemoved.contains(colObj.get(s)))) {
                sampleId++;
                for (int p = 0; p < dataset.rows(); ++p) {
                    newRawData[p][sampleId] = dataset.getMatrix().get(p, s);
                }
            }
        }

        for (String r : removeList) {
            dataset.hashCols.remove(r);
        }

        if ((dataset.rows() * newRawData[0].length) < (Integer.MAX_VALUE - 2)) {
            dataset = new SmallDoubleMatrixDataset<String, String>(new DenseDoubleMatrix2D(newRawData), dataset.hashRows, dataset.hashCols);
        } else {
            DenseLargeDoubleMatrix2D matrix = new DenseLargeDoubleMatrix2D(dataset.rows(), newRawData[0].length);
            matrix.assign(newRawData);
            dataset = new LargeDoubleMatrixDataset<String, String>(matrix, dataset.hashRows, dataset.hashCols);
        }
        return(dataset);
    }

    /**
     * Remove samples on forehand
     *
     * @param dataset DoubleMatrixDataset containing the matrix of interest
     * @param samplesToBeRemoved ArrayList<String> with identifiers of probes
     * that should be removed
     */
    public static void RemoveSamples(umcg.genetica.math.matrix2.DoubleMatrixDataset<String, String> dataset, HashSet<String> samplesToBeRemoved) {
        RemoveColumns(dataset, samplesToBeRemoved);
    }
    
    /**
     * Order cols to a index
     *
     * @param dataset DoubleMatrixDataset Expression matrix
     * @param hashRowsToInclude Ids of rowss to include
     */
    public static void ReorderCols(umcg.genetica.math.matrix2.DoubleMatrixDataset<String, String> dataset, LinkedHashMap<String, Integer> mappingIndex) {
        if(dataset instanceof SmallDoubleMatrixDataset){
            ReorderColsSmall((SmallDoubleMatrixDataset)dataset, mappingIndex);
        } else if(dataset instanceof LargeDoubleMatrixDataset){
            ReorderColsLarge((LargeDoubleMatrixDataset)dataset, mappingIndex);
        } else{
            throw new UnsupportedOperationException("Type of matrix not supported.");
        }
    }
    
    /**
     * Order cols to a index
     *
     * @param dataset DoubleMatrixDataset Expression matrix
     * @param hashRowsToInclude Ids of rowss to include
     */
    public static void ReorderColsSmall(umcg.genetica.math.matrix2.SmallDoubleMatrixDataset<String, String> dataset, LinkedHashMap<String, Integer> mappingIndex) {
        DenseDoubleMatrix2D newRawData = new DenseDoubleMatrix2D(dataset.rows(), dataset.columns());
        
        for (Entry<String, Integer> ent : mappingIndex.entrySet() ){
            int pos = dataset.getHashCols().get(ent.getKey());
            for (int p = 0; p < dataset.rows(); ++p) {
                newRawData.set(p, ent.getValue(), dataset.get(p, pos));
            }
        }
        dataset.setHashCols(mappingIndex);
        dataset.setMatrix(newRawData);
    }

    /**
     * Order cols to a index
     *
     * @param dataset DoubleMatrixDataset Expression matrix
     * @param hashRowsToInclude Ids of rowss to include
     */
    public static void ReorderColsLarge(umcg.genetica.math.matrix2.LargeDoubleMatrixDataset<String, String> dataset, LinkedHashMap<String, Integer> mappingIndex) {
        DenseLargeDoubleMatrix2D newRawData = new DenseLargeDoubleMatrix2D(dataset.rows(), dataset.columns());
        
        for (Entry<String, Integer> ent : mappingIndex.entrySet() ){
            int pos = dataset.getHashCols().get(ent.getKey());
            for (int p = 0; p < dataset.rows(); ++p) {
                newRawData.set(p, ent.getValue(), dataset.get(p, pos));
            }
        }
        
        dataset.setHashCols(mappingIndex);
        dataset.setMatrix(newRawData);
    }
    
    /**
     * Order columns
     *
     * @param dataset DoubleMatrixDataset Expression matrix
     */
    public static void OrderOnColumns(umcg.genetica.math.matrix2.DoubleMatrixDataset doubleMatrixDataset) {
        LinkedHashMap<String, Integer> newColHash = new LinkedHashMap<String, Integer>((int) Math.ceil(doubleMatrixDataset.columns() / 0.75));
        ArrayList<String> names = doubleMatrixDataset.getColObjects();
        Collections.sort(names);
        
        int pos = 0;
        for(String name : names){
            newColHash.put(name, pos);
            pos++;
        }
        ReorderCols(doubleMatrixDataset, newColHash);

    }

    /**
     * Order rows to a index
     *
     * @param dataset DoubleMatrixDataset Expression matrix
     * @param hashRowsToInclude Ids of rowss to include
     */
    public static void ReorderRows(umcg.genetica.math.matrix2.DoubleMatrixDataset<String, String> dataset, LinkedHashMap<String, Integer> mappingIndex) {
        if(dataset instanceof SmallDoubleMatrixDataset){
            ReorderRowsSmall( (SmallDoubleMatrixDataset) dataset, mappingIndex);
        } else if(dataset instanceof LargeDoubleMatrixDataset){
            ReorderRowsLarge((LargeDoubleMatrixDataset) dataset, mappingIndex);
        } else{
            throw new UnsupportedOperationException("Type of matrix not supported.");
        }
    }
    
    /**
     * Order rows to a index
     *
     * @param dataset DoubleMatrixDataset Expression matrix
     * @param hashRowsToInclude Ids of rowss to include
     */
    public static void ReorderRowsSmall(umcg.genetica.math.matrix2.SmallDoubleMatrixDataset<String, String> dataset, LinkedHashMap<String, Integer> mappingIndex) {
        DenseDoubleMatrix2D newRawData = new DenseDoubleMatrix2D(dataset.rows(), dataset.columns());
        
        for (Entry<String, Integer> ent : mappingIndex.entrySet() ){
            int pos = dataset.getHashRows().get(ent.getKey());
            for (int s = 0; s < dataset.columns(); ++s) {
                newRawData.set(ent.getValue(), s, dataset.get(pos, s));
            }
        }
        dataset.setHashRows(mappingIndex);
        dataset.setMatrix(newRawData);
    }

    /**
     * Order rows to a index
     *
     * @param dataset DoubleMatrixDataset Expression matrix
     * @param hashRowsToInclude Ids of rowss to include
     */
    public static void ReorderRowsLarge(umcg.genetica.math.matrix2.LargeDoubleMatrixDataset<String, String> dataset, LinkedHashMap<String, Integer> mappingIndex) {
        DenseLargeDoubleMatrix2D newRawData = new DenseLargeDoubleMatrix2D(dataset.rows(), dataset.columns());
        
        for (Entry<String, Integer> ent : mappingIndex.entrySet() ){
            int pos = dataset.getHashRows().get(ent.getKey());
            for (int s = 0; s < dataset.columns(); ++s) {
                newRawData.set(ent.getValue(), s, dataset.get(pos, s));
            }
        }
        
        dataset.setHashRows(mappingIndex);
        dataset.setMatrix(newRawData);
    }
    
    /**
     * Order rows
     *
     * @param dataset DoubleMatrixDataset Expression matrix
     */
    public static void OrderOnRows(umcg.genetica.math.matrix2.DoubleMatrixDataset doubleMatrixDataset) {
        
        LinkedHashMap<String, Integer> newRowHash = new LinkedHashMap<String, Integer>((int) Math.ceil(doubleMatrixDataset.rows() / 0.75));
        ArrayList<String> names = doubleMatrixDataset.getRowObjects();
        Collections.sort(names);
        
        int pos = 0;
        for(String name : names){
            newRowHash.put(name, pos);
            pos++;
        }
        ReorderRows(doubleMatrixDataset, newRowHash);

    }
    
    /**
     * Replace missing values in the double matrix per sample. Using either the
     * mean if useMedian is false or the median is useMedian is true.
     *
     * @param rawData
     * @param useMedian
     * @param NaValue
     */
    public static void ReplaceMissingValues(DoubleMatrix2D rawData, boolean useMedian, double NaValue) {

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
     * Remove rows without correct mapping known on forehand
     *
     * @param dataset DoubleMatrixDataset containing the matrix of interest
     * @param probesToKeep ArrayList<String> with identifiers of probes that
     * should be removed
     */
    public static DoubleMatrixDataset<String, String> FilterRows(umcg.genetica.math.matrix2.DoubleMatrixDataset<String, String> dataset, HashSet<String> probesToKeep) {
        int newSize = 0;
        HashSet<String> removeList = new HashSet<String>();

        for (String t : dataset.getRowObjects()) {
            if (probesToKeep.contains(t)) {
                newSize++;
            } else {
                removeList.add(t);
            }
        }

        double[][] newRawData = new double[newSize][dataset.columns()];

        int probeId = -1;
        ArrayList<String> rowObj = dataset.getRowObjects();

        for (int p = 0; p < dataset.rows(); ++p) {
            if ((probesToKeep.contains(rowObj.get(p)))) {
                probeId++;
                for (int s = 0; s < dataset.columns(); ++s) {
                    newRawData[probeId][s] = dataset.getMatrix().get(p, s);
                }
            }
        }

        for (String r : removeList) {
            dataset.hashRows.remove(r);
        }

        if ((dataset.columns() * newRawData.length) < (Integer.MAX_VALUE - 2)) {
            dataset = new SmallDoubleMatrixDataset<String, String>(new DenseDoubleMatrix2D(newRawData), dataset.hashRows, dataset.hashCols);
        } else {
            DenseLargeDoubleMatrix2D matrix = new DenseLargeDoubleMatrix2D(newRawData.length, dataset.columns());
            matrix.assign(newRawData);
            dataset = new LargeDoubleMatrixDataset<String, String>(matrix, dataset.hashRows, dataset.hashCols);
        }
        return(dataset);
    }

    /**
     * Filter out columns.
     * Keep all columns that in in the hashset.
     *
     * @param dataset
     * @param keepCols
     */
    public static DoubleMatrixDataset<String, String> FilterCols(umcg.genetica.math.matrix2.DoubleMatrixDataset<String, String> dataset, HashSet<String> keepCols) {
        int newSize = 0;
        HashSet<String> removeList = new HashSet<String>();

        for (String t : dataset.getColObjects()) {
            if (keepCols.contains(t)) {
                newSize++;
            } else {
                removeList.add(t);
            }
        }

        double[][] newRawData = new double[dataset.rows()][newSize];

        int sampleId = -1;

        ArrayList<String> colObj = dataset.getColObjects();

        for (int s = 0; s < dataset.columns(); ++s) {
            if ((keepCols.contains(colObj.get(s)))) {
                sampleId++;
                for (int p = 0; p < dataset.rows(); ++p) {
                    newRawData[p][sampleId] = dataset.getMatrix().get(p, s);
                }
            }
        }

        for (String r : removeList) {
            dataset.hashCols.remove(r);
        }

        if ((dataset.rows() * newRawData[0].length) < (Integer.MAX_VALUE - 2)) {
            dataset = new SmallDoubleMatrixDataset<String, String>(new DenseDoubleMatrix2D(newRawData), dataset.hashRows, dataset.hashCols);
        } else {
            DenseLargeDoubleMatrix2D matrix = new DenseLargeDoubleMatrix2D(dataset.rows(), newRawData[0].length);
            matrix.assign(newRawData);
            dataset = new LargeDoubleMatrixDataset<String, String>(matrix, dataset.hashRows, dataset.hashCols);
        }
        return(dataset);
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
}
