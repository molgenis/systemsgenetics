/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.io;

import java.io.File;
import java.io.FileFilter;
import java.io.IOException;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;
import org.apache.commons.io.filefilter.WildcardFileFilter;
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

		if (!blockReader.getColumnIdentifiers().containsAll(identifiers)) {
			throw new IOException("Not all requested genes found in columns");
		}

		if (!blockReader.getRowIdentifiers().containsAll(identifiers)) {
			throw new IOException("Not all requested genes found in rows");
		}

		return blockReader.loadSubsetOfRows(identifiers).viewColSelection(identifiers);

	}

	@Override
	public Set<String> getGenes() throws IOException {

		File x = new File(prefix);
		
		File dir = x.getParentFile();
		FileFilter fileFilter = new WildcardFileFilter(x.getName() + "*" + suffix + ".datg");
		File[] files = dir.listFiles(fileFilter);
		
		HashSet<String> allGenes = new HashSet<>();
		
		for(File file : files){
			allGenes.addAll(new DoubleMatrixDatasetRowCompressedReader(file.getAbsolutePath()).getRowIdentifiers());
		}

		return allGenes;
	}

}
