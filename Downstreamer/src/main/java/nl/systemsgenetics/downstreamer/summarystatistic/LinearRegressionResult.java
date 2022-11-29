package nl.systemsgenetics.downstreamer.summarystatistic;

import umcg.genetica.math.matrix2.DoubleMatrixDataset;

public class LinearRegressionResult {

    private final DoubleMatrixDataset<String, String> beta;
    private final DoubleMatrixDataset<String, String> standardError;
    private final DoubleMatrixDataset<String, String> explainedVariance;
    private final int degreesOfFreedom;

    public LinearRegressionResult(DoubleMatrixDataset<String, String> beta, DoubleMatrixDataset<String, String> standardError, int degreesOfFreedom) {
        this.beta = beta;
        this.standardError = standardError;
        this.degreesOfFreedom = degreesOfFreedom;
        this.explainedVariance = null;
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
}
