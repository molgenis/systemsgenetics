/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.io;

import java.io.IOException;
import java.util.Collection;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRowCompressedReader;

/**
 *
 * @author patri
 */
public class BlockPerFileDiagonalDoubleMatrixProvider implements BlockDiagonalDoubleMatrixProvider {

	final String prefix;
	final String suffix;

	public BlockPerFileDiagonalDoubleMatrixProvider(String prefix, String suffix) {
		this.prefix = prefix;
		this.suffix = suffix;
	}

	@Override
	public DoubleMatrixDataset<String, String> viewBlock(String block, Collection<String> identifiers) throws IOException, Exception {

		DoubleMatrixDatasetRowCompressedReader blockReader = new DoubleMatrixDatasetRowCompressedReader(prefix + block + suffix);
		
		if(!blockReader.getColumnIdentifiers().containsAll(identifiers)){
			throw new IOException("Not all requested genes found in columns");
		}
		
		if(!blockReader.getRowIdentifiers().containsAll(identifiers)){
			throw new IOException("Not all requested genes found in rows");
		}
		
		return blockReader.loadSubsetOfRows(identifiers).viewColSelection(identifiers);

	}

}
