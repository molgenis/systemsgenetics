package nl.systemsgenetics.downstreamer.io;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.util.List;

public interface BlockDiagonalDoubleMatrixProvider {

    DoubleMatrix2D viewBlock(int[] indices);
    DoubleMatrixDataset<String, String> viewBlock(List<String> items);
    DoubleMatrixDataset<String, String> getAsFullMatrix();

    List<String> getRowNames();
    List<String> getCollumnNames();

    int rows();
    int columns();

}
