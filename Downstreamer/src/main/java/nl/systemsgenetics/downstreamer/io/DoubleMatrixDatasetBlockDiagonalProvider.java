package nl.systemsgenetics.downstreamer.io;

import java.util.Collection;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.util.List;

public class DoubleMatrixDatasetBlockDiagonalProvider implements  BlockDiagonalDoubleMatrixProvider {

    private final DoubleMatrixDataset<String, String> data;

    public DoubleMatrixDatasetBlockDiagonalProvider(DoubleMatrixDataset<String, String> data) {
        this.data = data;
    }

    @Override
    public DoubleMatrixDataset<String, String> viewBlock(String block, Collection<String> items) {
        return data.viewSelection(items, items);
    }

    public List<String> getRowNames() {
        return data.getRowObjects();
    }

    public List<String> getCollumnNames() {
        return data.getColObjects();
    }

    public int rows() {
        return data.rows();
    }

    public int columns() {
        return data.columns();
    }


}
