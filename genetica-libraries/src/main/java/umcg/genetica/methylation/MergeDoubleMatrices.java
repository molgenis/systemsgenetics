/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.methylation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author MarcJan
 */
public class MergeDoubleMatrices {

    public static DoubleMatrixDataset<String, String> combine(ArrayList<DoubleMatrixDataset> datasets) {

        HashSet<String> removeList = new HashSet<String>();

        for (int i = 0; i < (datasets.size()); ++i) {
            for (int j = 0; j < (datasets.size()); ++j) {
                if (i != j) {
                    DoubleMatrixDataset<String, String> ds1 = datasets.get(i);
                    DoubleMatrixDataset<String, String> ds2 = datasets.get(j);
                    for (Entry<String, Integer> t : ds1.hashRows.entrySet()) {
                        if (!(ds2.hashRows.containsKey(t.getKey()))) {
                            removeList.add(t.getKey());
                        }
                    }
                }
            }
        }

        System.out.println("Number of probes missing in one or more datasets: " + removeList.size());

        //ToDo handle duplicate entry names. 

        HashMap<String, Integer> sampleToDataset = new HashMap<String, Integer>();
        ArrayList<String> newColumnNames = new ArrayList<String>();

        for (int i = 0; i < (datasets.size()); ++i) {
            if (removeList.size() > 0) {
                CleanDoubleMatrix.RemoveRows(datasets.get(i), removeList);
            }
            
            CleanDoubleMatrix.OrderOnRows(datasets.get(i));

            HashMap<String, Integer> t = (HashMap<String, Integer>) datasets.get(i).hashCols;
            for (Entry<String, Integer> cols : t.entrySet()) {
                if (!sampleToDataset.containsKey(cols.getKey())) {
                    sampleToDataset.put(cols.getKey(), i);
                    newColumnNames.add(cols.getKey());
                }
            }
        }

        System.out.println("Number of unique columnnames : " + sampleToDataset.size());
        System.out.println("Number of columnnames : " + newColumnNames.size());

        double[][] newRawData = new double[datasets.get(0).nrRows][newColumnNames.size()];

        for (int p = 0; p < datasets.get(0).nrRows; ++p) {
            for (int s = 0; s < newColumnNames.size(); ++s) {
                int dataSetNr = sampleToDataset.get(newColumnNames.get(s));

                //We shouldnt have to do this conversion.
                int originalS = (Integer) datasets.get(dataSetNr).hashCols.get(newColumnNames.get(s));
                newRawData[p][s] = datasets.get(dataSetNr).rawData[p][originalS];
            }
        }

        DoubleMatrixDataset<String, String> newDataset = new DoubleMatrixDataset<String, String>();
        newDataset.rowObjects = datasets.get(0).rowObjects;
        newDataset.nrRows = datasets.get(0).nrRows;

        newDataset.colObjects = newColumnNames;
        newDataset.nrCols = sampleToDataset.size();

        newDataset.rawData = newRawData;
        newDataset.recalculateHashMaps();

        return (newDataset);
    }
}
