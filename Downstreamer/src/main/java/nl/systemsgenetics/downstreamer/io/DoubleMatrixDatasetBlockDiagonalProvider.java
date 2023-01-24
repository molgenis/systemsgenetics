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
