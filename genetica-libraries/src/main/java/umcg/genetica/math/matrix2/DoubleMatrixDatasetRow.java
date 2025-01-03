package umcg.genetica.math.matrix2;

public class DoubleMatrixDatasetRow {

    private final double[] data;
    private final String rowName;

    public DoubleMatrixDatasetRow(double[] data, String rowName) {
        this.data = data;
        this.rowName = rowName;
    }

    public double[] getData() {
        return data;
    }

    public String getRowName() {
        return rowName;
    }

}
