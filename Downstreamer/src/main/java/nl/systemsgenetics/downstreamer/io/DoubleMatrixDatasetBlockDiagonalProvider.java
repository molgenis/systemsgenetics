package nl.systemsgenetics.downstreamer.io;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import org.jblas.DoubleMatrix;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.util.List;

public class DoubleMatrixDatasetBlockDiagonalProvider implements  BlockDiagonalDoubleMatrixProvider {

    private final DoubleMatrixDataset<String, String> data;

    public DoubleMatrixDatasetBlockDiagonalProvider(DoubleMatrixDataset<String, String> data) {
        this.data = data;
    }

    @Override
    public DoubleMatrix2D viewBlock(int[] indices) {
        return  data.getMatrix().viewSelection(indices, indices);
    }

    @Override
    public DoubleMatrixDataset<String, String> viewBlock(List<String> items) {
        return data.viewSelection(items, items);
    }

    @Override
    public DoubleMatrixDataset<String, String> getAsFullMatrix() {
        return data;
    }

    @Override
    public List<String> getRowNames() {
        return data.getRowObjects();
    }

    @Override
    public List<String> getCollumnNames() {
        return data.getColObjects();
    }

    @Override
    public int rows() {
        return data.rows();
    }

    @Override
    public int columns() {
        return data.columns();
    }


}
