/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.matrix2;

import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;

/**
 *
 * @author patri
 */
public class DoubleMatrixDatasetFastSubsetLoader {

	private final LinkedHashMap<String, Integer> originalRowMap;
	private final LinkedHashMap<String, Integer> originalColMap;
	private final String fileName;

	public DoubleMatrixDatasetFastSubsetLoader(String fileName) throws IOException {
		this.fileName = fileName;

		originalRowMap = DoubleMatrixDataset.loadIdentifiers(fileName + ".rows.txt");
		originalColMap = DoubleMatrixDataset.loadIdentifiers(fileName + ".cols.txt");

	}

	/**
	 *
	 * @param rowsToView
	 * @return subset of rows in order of rowsToView
	 * @throws IOException
	 */
	public DoubleMatrixDataset<String, String> loadSubsetOfRowsBinaryDoubleData(LinkedHashSet<String> rowsToView) throws IOException {
		return DoubleMatrixDataset.loadSubsetOfRowsBinaryDoubleData(fileName, rowsToView, originalRowMap, originalColMap);
	}
	
	/**
	 *
	 * @param rowsToView
	 * @return subset of rows in order of rowsToView
	 * @throws IOException
	 */
	public DoubleMatrixDataset<String, String> loadSubsetOfRowsBinaryDoubleData(String[] rowsToView) throws IOException {

		LinkedHashSet<String> rowsToViewHash = new LinkedHashSet<>(rowsToView.length);

		for (String rowToView : rowsToView) {
			rowsToViewHash.add(rowToView);
		}

		if (rowsToViewHash.size() != rowsToView.length) {
			throw new RuntimeException("Duplicates in rows not allowed");
		}

		return loadSubsetOfRowsBinaryDoubleData(rowsToViewHash);

	}

}
