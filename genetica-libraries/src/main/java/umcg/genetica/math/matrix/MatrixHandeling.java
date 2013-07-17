/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.matrix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
import org.apache.commons.collections.primitives.ArrayDoubleList;

/**
 *
 * @author Marc Jan
 */
public class MatrixHandeling {

    /**
     * Remove columns with to many missing values
     *
     * @param dataset DoubleMatrixDataset Expression matrix
     * @param hashColumnsToInclude Ids of samples to include
     */
    public static void RemoveColumnsWithToManyMissingValues(DoubleMatrixDataset<String, String> dataset, int maxMissingValuesPerColumn) {

        ArrayList<String> columnsToInclude = new ArrayList<String>();

        for (int c = 0; c < dataset.nrCols; ++c) {
            int nrMissing = 0;
            for (int r = 0; r < dataset.nrRows; ++r) {
                if (dataset.rawData[r][c] == -999) {
                    nrMissing++;
                }
            }
            if (nrMissing >= maxMissingValuesPerColumn) {
                System.out.println("Excluding:\t" + c + "\t" + dataset.colObjects.get(c) + "\t" + nrMissing);
            } else if (nrMissing > 0) {
                //System.out.println("Missing values in:\t" + c + "\t" + dataset.colObjects.get(c) + "\t" + nrMissing);
                columnsToInclude.add(dataset.colObjects.get(c));
            } else {
                columnsToInclude.add(dataset.colObjects.get(c));
            }
        }

        System.out.println(columnsToInclude.size() + "\t" + dataset.nrCols + "\t" + dataset.colObjects.size());

        double[][] newRawData = new double[dataset.nrRows][columnsToInclude.size()];
        String[] newColumnNames = new String[columnsToInclude.size()];

        int sampleId = -1;

        for (int s = 0; s < dataset.nrCols; ++s) {
            if (columnsToInclude.contains(dataset.colObjects.get(s))) {
                sampleId++;
                newColumnNames[sampleId] = dataset.colObjects.get(s);
                for (int p = 0; p < dataset.nrRows; ++p) {
                    newRawData[p][sampleId] = dataset.rawData[p][s];
                }
            }
        }

        dataset.colObjects = Arrays.asList(newColumnNames);
        dataset.nrCols = columnsToInclude.size();
        dataset.rawData = newRawData;

    }

    /**
     * Remove rows with to many missing values
     *
     * @param dataset DoubleMatrixDataset Expression matrix
     * @param hashRowsToInclude Ids of rowss to include
     */
    public static void RemoveRowsWithToManyMissingValues(DoubleMatrixDataset<String, String> dataset, int maxMissingValuesPerRow) {

        String[] rowNames = dataset.rowObjects.toArray(new String[0]);

        HashMap<String, Boolean> hashRowsToInclude = new HashMap<String, Boolean>();

        for (int r = 0; r < dataset.nrRows; ++r) {
            int nrMissing = 0;
            for (int c = 0; c < dataset.nrCols; ++c) {
                if (dataset.rawData[r][c] == -999) {
                    nrMissing++;
                }
            }
            if (nrMissing >= maxMissingValuesPerRow) {
                System.out.println("Excluding:\t" + r + "\t" + rowNames[r] + "\t" + nrMissing);
            } else if (nrMissing > 0) {
                //System.out.println("Missing values in:\t" + r + "\t" + rowNames[r] + "\t" + nrMissing);
                hashRowsToInclude.put(rowNames[r], null);
            } else {
                hashRowsToInclude.put(rowNames[r], null);
            }
        }

        double[][] newRawData = new double[hashRowsToInclude.size()][dataset.nrCols];
        String[] newRowNames = new String[hashRowsToInclude.size()];

        int probeId = -1;

        for (int p = 0; p < dataset.nrRows; ++p) {
            if (hashRowsToInclude.containsKey(rowNames[p])) {
                probeId++;
                newRowNames[probeId] = rowNames[p];
                for (int s = 0; s < dataset.nrCols; ++s) {
                    newRawData[probeId][s] = dataset.rawData[p][s];
                }
            }
        }

        dataset.rowObjects = Arrays.asList(newRowNames);
        dataset.nrRows = hashRowsToInclude.size();
        dataset.rawData = newRawData;
        dataset.recalculateHashMaps();
    }

    /**
     * Remove identical samples, based on name and/or based on all expression
     * value (advanced = true)
     *
     * @param dataset DoubleMatrixDataset Expression matrix
     * @param advanced Do sample check based on double values
     */
    public static void RemoveDuplicatesSamples(DoubleMatrixDataset<String, String> dataset, boolean advanced) {

        if (advanced) {
            throw new UnsupportedOperationException("Not yet implemented");
        }

        ArrayList<Integer> removeEntry = new ArrayList<Integer>();

        for (int i = 0; i < dataset.nrCols; ++i) {
            if (dataset.colObjects.get(i) == null) {
                removeEntry.add(i);
            } else {
                List<String> tmp = dataset.colObjects.subList(i + 1, dataset.nrCols);
                if (tmp.contains(dataset.colObjects.get(i))) {
                    removeEntry.add(i);
                }
            }
        }

        if (removeEntry.size() > 0) {
            int newSize = (dataset.nrCols - removeEntry.size());
            double[][] newRawData = new double[dataset.nrRows][newSize];
            String[] newColumnNames = new String[newSize];

            int sampleId = -1;

            for (int s = 0; s < dataset.nrCols; ++s) {
                if (!removeEntry.contains(s)) {
                    sampleId++;
                    newColumnNames[sampleId] = dataset.colObjects.get(s);
                    for (int p = 0; p < dataset.nrRows; ++p) {
                        newRawData[p][sampleId] = dataset.rawData[p][s];
                    }
                }
            }

            dataset.colObjects = Arrays.asList(newColumnNames);
            dataset.nrCols = newSize;
            dataset.rawData = newRawData;
            dataset.recalculateHashMaps();
        }
    }

    /**
     * Remove rows without correct mapping known on forehand
     *
     * @param dataset DoubleMatrixDataset containing the matrix of interest
     * @param probesToBeRemoved ArrayList<String> with identifiers of probes
     * that should be removed
     */
    public static void RemoveRows(DoubleMatrixDataset<String, String> dataset, HashSet<String> probesToBeRemoved) {
        int newSize = 0;
        for (String t : dataset.rowObjects) {
            if (!probesToBeRemoved.contains(t)) {
                newSize++;
            }
        }


        double[][] newRawData = new double[newSize][dataset.nrCols];
        String[] newRowNames = new String[newSize];

        int probeId = -1;

        for (int p = 0; p < dataset.nrRows; ++p) {
            if (!(probesToBeRemoved.contains(dataset.rowObjects.get(p)))) {
                probeId++;
                newRowNames[probeId] = dataset.rowObjects.get(p);

                newRawData[probeId] = dataset.rawData[p];

            }
        }

        dataset.rowObjects = Arrays.asList(newRowNames);
        dataset.nrRows = newSize;
        dataset.rawData = newRawData;
        dataset.recalculateHashMaps();
    }

    /**
     * Remove rows without correct mapping known on forehand
     *
     * @param dataset DoubleMatrixDataset containing the matrix of interest
     * @param probesToKeep ArrayList<String> with identifiers of probes that
     * should be removed
     */
    public static void FilterRows(DoubleMatrixDataset<String, String> dataset, HashSet<String> probesToKeep) {
        int newSize = 0;
        for (String t : dataset.rowObjects) {
            if (probesToKeep.contains(t)) {
                newSize++;
            }
        }

        double[][] newRawData = new double[newSize][dataset.nrCols];
        String[] newRowNames = new String[newSize];

        int probeId = -1;

        for (int p = 0; p < dataset.nrRows; ++p) {
            if ((probesToKeep.contains(dataset.rowObjects.get(p)))) {
                probeId++;
                newRowNames[probeId] = dataset.rowObjects.get(p);
                for (int s = 0; s < dataset.nrCols; ++s) {
                    newRawData[probeId][s] = dataset.rawData[p][s];
                }
            }
        }

        dataset.rowObjects = Arrays.asList(newRowNames);
        dataset.nrRows = newSize;
        dataset.rawData = newRawData;
        dataset.recalculateHashMaps();
    }

    /**
     * Remove probes without correct mapping known on forehand
     *
     * @param dataset DoubleMatrixDataset containing the matrix of interest
     * @param probesToBeRemoved ArrayList<String> with identifiers of probes
     * that should be removed
     */
    public static void RemoveProbes(DoubleMatrixDataset<String, String> dataset, HashSet<String> probesToBeRemoved) {
        RemoveRows(dataset, probesToBeRemoved);
    }

    /**
     * Remove samples on forehand
     *
     * @param dataset DoubleMatrixDataset containing the matrix of interest
     * @param samplesToBeRemoved ArrayList<String> with identifiers of probes
     * that should be removed
     */
    public static void RemoveColumns(DoubleMatrixDataset<String, String> dataset, HashSet<String> samplesToBeRemoved) {
        int newSize = 0;
        for (String t : dataset.colObjects) {
            if (!samplesToBeRemoved.contains(t)) {
                newSize++;
            }
        }

        double[][] newRawData = new double[dataset.nrRows][newSize];
        String[] newColNames = new String[newSize];

        int sampleId = -1;

        for (int s = 0; s < dataset.nrCols; ++s) {
            if (!(samplesToBeRemoved.contains(dataset.colObjects.get(s)))) {
                sampleId++;
                newColNames[sampleId] = dataset.colObjects.get(s);
                for (int p = 0; p < dataset.nrRows; ++p) {
                    newRawData[p][sampleId] = dataset.rawData[p][s];
                }
            }
        }

        dataset.colObjects = Arrays.asList(newColNames);
        dataset.nrCols = newSize;
        dataset.rawData = newRawData;
        dataset.recalculateHashMaps();
    }

    /**
     * Remove samples on forehand
     *
     * @param dataset DoubleMatrixDataset containing the matrix of interest
     * @param samplesToBeRemoved ArrayList<String> with identifiers of probes
     * that should be removed
     */
    public static void RemoveSamples(DoubleMatrixDataset<String, String> dataset, HashSet<String> samplesToBeRemoved) {
        RemoveColumns(dataset, samplesToBeRemoved);
    }

    /**
     * Order rows
     *
     * @param dataset DoubleMatrixDataset Expression matrix
     */
    public static void OrderOnRows(DoubleMatrixDataset doubleMatrixDataset) {
        HashMap<String, Integer> originalRowMap = (HashMap<String, Integer>) doubleMatrixDataset.hashRows;
        ArrayList<String> rowNames = new ArrayList<String>(doubleMatrixDataset.rowObjects);
        HashMap<String, Integer> newRowMap = new HashMap<String, Integer>();

        Collections.sort(rowNames);

        for (int i = 0; i < rowNames.size(); ++i) {
            newRowMap.put(rowNames.get(i), i);
        }

        double[][] newRawData = new double[doubleMatrixDataset.nrRows][doubleMatrixDataset.nrCols];

        for (String rName : rowNames) {
            int originalprobeId = originalRowMap.get(rName);
            int newprobeId = newRowMap.get(rName);

            newRawData[newprobeId] = doubleMatrixDataset.rawData[originalprobeId];

        }

        doubleMatrixDataset.rowObjects = rowNames;
        doubleMatrixDataset.rawData = newRawData;
        doubleMatrixDataset.recalculateHashMaps();

    }

    /**
     * Order columns
     *
     * @param dataset DoubleMatrixDataset Expression matrix
     */
    public static void OrderOnColumns(DoubleMatrixDataset doubleMatrixDataset) {
        HashMap<String, Integer> originalColMap = (HashMap<String, Integer>) doubleMatrixDataset.hashCols;
        ArrayList<String> colNames = new ArrayList<String>(doubleMatrixDataset.colObjects);
        HashMap<String, Integer> newColMap = new HashMap<String, Integer>();

        Collections.sort(colNames);

        for (int i = 0; i < colNames.size(); ++i) {
            newColMap.put(colNames.get(i), i);
        }

        double[][] newRawData = new double[doubleMatrixDataset.nrRows][doubleMatrixDataset.nrCols];

        for (String cName : colNames) {
            int originalSampleId = originalColMap.get(cName);
            int newSampleId = newColMap.get(cName);

            for (int r = 0; r < doubleMatrixDataset.nrRows; ++r) {
                newRawData[r][newSampleId] = doubleMatrixDataset.rawData[r][originalSampleId];
            }
        }

        doubleMatrixDataset.colObjects = colNames;
        doubleMatrixDataset.rawData = newRawData;
        doubleMatrixDataset.recalculateHashMaps();

    }

    /**
     * Order rows to a index
     *
     * @param dataset DoubleMatrixDataset Expression matrix
     * @param hashRowsToInclude Ids of rowss to include
     */
    public static void ReorderRows(DoubleMatrixDataset<String, String> dataset, HashMap<String, Integer> mappingIndex) {

        double[][] newRawData = new double[mappingIndex.size()][dataset.nrCols];
        String[] newRowNames = new String[mappingIndex.size()];

        for (Entry<String, Integer> instance : mappingIndex.entrySet()) {
            newRowNames[instance.getValue()] = instance.getKey();
            for (int s = 0; s < dataset.nrCols; ++s) {
                newRawData[instance.getValue()][s] = dataset.rawData[dataset.hashRows.get(instance.getKey())][s];
            }
        }
        dataset.rowObjects = Arrays.asList(newRowNames);
        dataset.nrRows = mappingIndex.size();
        dataset.rawData = newRawData;
        dataset.recalculateHashMaps();
    }

    /**
     * Replace missing values in the double matrix per sample. Using either the
     * mean if useMedian is false or the median is useMedian is true.
     *
     * @param rawData
     * @param useMedian
     */
    public static void ReplaceMissingValues(double[][] rawData, boolean useMedian) {

        for (int s = 0; s < rawData[1].length; ++s) {

            System.out.println("Processing sample: " + s);
            boolean needsReplacement = false;
            ArrayDoubleList nonNAvalues = new ArrayDoubleList();

            for (int p = 0; p < rawData.length; ++p) {
                if (rawData[p][s] == -999) {
                    needsReplacement = true;
                } else {
                    nonNAvalues.add(rawData[p][s]);
                }
            }

            if (needsReplacement) {
                double replacementValue;
                if (useMedian) {
                    replacementValue = JSci.maths.ArrayMath.median(nonNAvalues.toArray(new double[0]));
                } else {
                    replacementValue = JSci.maths.ArrayMath.mean(nonNAvalues.toArray(new double[0]));

                }

                for (int p = 0; p < rawData.length; ++p) {
                    if (rawData[p][s] == -999) {
                        rawData[p][s] = replacementValue;
                    }
                }
            }
        }
    }

    /**
     * Method to apply all QC steps.
     *
     *
     * @param dataset
     * @param probesToBeRemoved
     * @param samplesToBeRemoved
     * @param replaceKnownOutOfRangeValues
     * @param maxMissingSamplesMissingAProbe
     * @param maxMissingProbesMissingPerSample
     */
    public static void performQC(DoubleMatrixDataset<String, String> dataset, ArrayList<String> probesToBeRemoved, ArrayList<String> samplesToBeRemoved, boolean replaceKnownOutOfRangeValues, int maxMissingSamplesMissingAProbe, int maxMissingProbesMissingPerSample) {

        int originalNumberSamples = dataset.nrCols;
        int originalNumberProbes = dataset.nrRows;

        //Remove probes withoutmapping to HG19

        RemoveDuplicatesSamples(dataset, false);
        System.out.println(dataset.nrCols);

        System.out.println(originalNumberSamples);
        RemoveProbes(dataset, new HashSet(probesToBeRemoved));
        System.out.println(originalNumberProbes - dataset.nrRows);

        RemoveSamples(dataset, new HashSet(probesToBeRemoved));

        ArrayList<String> extraSamplesToBeRemoved = checkMinAndMaxPerSample(dataset, replaceKnownOutOfRangeValues);

        RemoveSamples(dataset, new HashSet(probesToBeRemoved));

        RemoveColumnsWithToManyMissingValues(dataset, maxMissingProbesMissingPerSample);
        System.out.println("Number of samples removed: " + (originalNumberSamples - dataset.rawData[1].length));

        RemoveRowsWithToManyMissingValues(dataset, maxMissingSamplesMissingAProbe);
        System.out.println("Number of probes removed: " + (originalNumberProbes - dataset.rawData.length));


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
                        } else if (dataset.rawData[r][c] == (-3.4d * Math.pow(10, 38))) {
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

    public static void FilterCols(DoubleMatrixDataset<String, String> dataset, HashSet<String> keepCols) {
        int newSize = 0;
        for (String t : dataset.colObjects) {
            if (keepCols.contains(t)) {
                newSize++;
            }
        }

        double[][] newRawData = new double[dataset.nrRows][newSize];
        String[] newColNames = new String[newSize];

        int sampleId = -1;

        for (int s = 0; s < dataset.nrCols; ++s) {
            if ((keepCols.contains(dataset.colObjects.get(s)))) {
                sampleId++;
                newColNames[sampleId] = dataset.colObjects.get(s);
                for (int p = 0; p < dataset.nrRows; ++p) {
                    newRawData[p][sampleId] = dataset.rawData[p][s];
                }
            }
        }

        dataset.colObjects = Arrays.asList(newColNames);
        dataset.nrCols = newSize;
        dataset.rawData = newRawData;
        dataset.recalculateHashMaps();
    }
}
