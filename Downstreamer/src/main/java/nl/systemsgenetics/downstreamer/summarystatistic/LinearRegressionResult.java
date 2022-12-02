package nl.systemsgenetics.downstreamer.summarystatistic;

import org.apache.commons.lang3.NotImplementedException;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.util.Collection;

public class LinearRegressionResult {

    private final DoubleMatrixDataset<String, String> beta;
    private final DoubleMatrixDataset<String, String> standardError;
    private final DoubleMatrixDataset<String, String> explainedVariance;
    private final int degreesOfFreedom;

    private final DoubleMatrixDataset<String, String> tStatisticCache;
    private final DoubleMatrixDataset<String, String> pValueCache;

    public LinearRegressionResult(Collection<String> rownames, Collection<String> colnames, int degreesOfFreedom) {
        this.beta = new DoubleMatrixDataset<>(rownames, colnames);
        this.standardError = new DoubleMatrixDataset<>(rownames, colnames);
        this.explainedVariance = null;
        this.degreesOfFreedom = degreesOfFreedom;

        this.tStatisticCache = null;
        this.pValueCache = null;
    }

    public LinearRegressionResult(DoubleMatrixDataset<String, String> beta, DoubleMatrixDataset<String, String> standardError, int degreesOfFreedom) {
        this.beta = beta;
        this.standardError = standardError;
        this.degreesOfFreedom = degreesOfFreedom;
        this.explainedVariance = null;
        this.tStatisticCache = null;
        this.pValueCache = null;
    }

    public DoubleMatrixDataset<String, String> getBeta() {
        return beta;
    }

    public DoubleMatrixDataset<String, String> getStandardError() {
        return standardError;
    }

    public DoubleMatrixDataset<String, String> getExplainedVariance() {
        return explainedVariance;
    }

    public int getDegreesOfFreedom() {
        return degreesOfFreedom;
    }

    public void appendBetas(int rowNumber, double[] betas) {
        //int rowNumber = beta.getRowIndex(row);
        for (int i=0; i < betas.length; i++) {
            beta.setElementQuick(rowNumber, i, betas[i]);
        }
    }

    public void appendSe(int rowNumber, double[] se) {
        //int rowNumber = standardError.getRowIndex(row);
        for (int i=0; i<se.length; i++) {
            standardError.setElementQuick(rowNumber, i, se[i]);
        }
    }

    public DoubleMatrixDataset<String, String> getTstats() {
        throw new NotImplementedException("Not yet implemented");
    }

    public DoubleMatrixDataset<String, String> getPvalues() {
        throw new NotImplementedException("Not yet implemented");
    }


    public void save(String basePath, boolean isBinary) throws Exception{

        String finalPath = basePath + "_df_" + degreesOfFreedom;

        if (isBinary) {
            beta.saveBinary(finalPath + "_betas");
            standardError.saveBinary(finalPath + "_se");
        } else {
            beta.save(finalPath + "_betas");
            standardError.save(finalPath + "_se");
        }
    }
}
