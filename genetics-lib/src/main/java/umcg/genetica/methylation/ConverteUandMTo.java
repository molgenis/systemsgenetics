/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.methylation;

import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author MarcJan
 */
public class ConverteUandMTo {

    public static DoubleMatrixDataset<String, String> Mvalue(DoubleMatrixDataset<String, String> uChanel, DoubleMatrixDataset<String, String> mChanel) {
        DoubleMatrixDataset<String, String> mVals = new DoubleMatrixDataset<String, String>();

        double multiplier = 1.0d / Math.log10(2.0d);

        double[][] uRaw = uChanel.rawData;
        double[][] mRaw = mChanel.rawData;

        int probeCount = uRaw.length;
        int sampleCount = uRaw[probeCount - 1].length;

        double[][] raw = new double[probeCount][sampleCount];

        for (int p = 0; p < probeCount; p++) {
            for (int s = 0; s < sampleCount; s++) {
                if (uRaw[p][s] < 0) {
                    uRaw[p][s] = 0;
                }
                if (mRaw[p][s] < 0) {
                    mRaw[p][s] = 0;
                }
                raw[p][s] = Math.log10(((mRaw[p][s] + 100.0d) / (uRaw[p][s] + 100.0d))) * multiplier;
            }
        }

        mVals.rawData = raw;
        mVals.colObjects = uChanel.colObjects;
        mVals.rowObjects = uChanel.rowObjects;
        mVals.nrCols = uChanel.nrCols;
        mVals.nrRows = uChanel.nrRows;
        mVals.recalculateHashMaps();

        return (mVals);
    }

    public static DoubleMatrixDataset<String, String> Beta(DoubleMatrixDataset<String, String> uChanel, DoubleMatrixDataset<String, String> mChanel) {
        DoubleMatrixDataset<String, String> beta = new DoubleMatrixDataset<String, String>();

        double[][] uRaw = uChanel.rawData;
        double[][] mRaw = mChanel.rawData;

        int probeCount = uRaw.length;
        int sampleCount = uRaw[probeCount - 1].length;

        double[][] raw = new double[probeCount][sampleCount];

        for (int p = 0; p < probeCount; p++) {
            for (int s = 0; s < sampleCount; s++) {
                if (uRaw[p][s] < 0) {
                    uRaw[p][s] = 0;
                }
                if (mRaw[p][s] < 0) {
                    mRaw[p][s] = 0;
                }
                raw[p][s] = (mRaw[p][s] / (mRaw[p][s] + uRaw[p][s] + 100.0d));
            }
        }

        beta.rawData = raw;
        beta.colObjects = uChanel.colObjects;
        beta.rowObjects = uChanel.rowObjects;
        beta.recalculateHashMaps();

        return (beta);
    }
    
}
