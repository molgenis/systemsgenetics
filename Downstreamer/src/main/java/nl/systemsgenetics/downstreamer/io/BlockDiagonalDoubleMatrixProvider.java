package nl.systemsgenetics.downstreamer.io;

import java.io.IOException;
import java.util.Collection;
import java.util.Set;

import umcg.genetica.math.matrix2.DoubleMatrixDataset;

public interface BlockDiagonalDoubleMatrixProvider {

    DoubleMatrixDataset<String, String> viewBlock(String block, Collection<String> identifiers) throws IOException, Exception;

    Set<String> getGenes() throws IOException;
//    List<String> getCollumnNames();
//
//    int rows();
//    int columns();

}
