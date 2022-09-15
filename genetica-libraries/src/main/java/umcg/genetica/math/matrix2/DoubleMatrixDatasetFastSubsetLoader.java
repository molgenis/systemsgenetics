/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.matrix2;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map;

/**
 * @author patri
 */
public class DoubleMatrixDatasetFastSubsetLoader {

	private LinkedHashMap<String, Integer> originalRowMap;
	private LinkedHashMap<String, Integer> originalColMap;
	private final String fileName;

	public DoubleMatrixDatasetFastSubsetLoader(String fileName) throws IOException {
		this.fileName = fileName;

		// strip .dat or .dat.gz from filename
		if (fileName.endsWith(".dat")) {
			fileName = fileName.substring(0, fileName.length() - 4);
		} else if (fileName.endsWith(".dat.gz")) {
			fileName = fileName.substring(0, fileName.length() - 7);
		}

		//Now load the row and column identifiers from files
		originalRowMap = null;
		originalColMap = null;
		if (new File(fileName + ".rows.txt").exists()) {
			originalRowMap = DoubleMatrixDataset.loadIdentifiers(fileName + ".rows.txt");
		} else if (new File(fileName + ".rows.txt.gz").exists()) {
			originalRowMap = DoubleMatrixDataset.loadIdentifiers(fileName + ".rows.txt.gz");
		} else {
			throw new FileNotFoundException("File not found: " + fileName + ".rows.txt or " + fileName + ".rows.txt.gz");
		}

		if (new File(fileName + ".cols.txt").exists()) {
			originalColMap = DoubleMatrixDataset.loadIdentifiers(fileName + ".cols.txt");
		} else if (new File(fileName + ".cols.txt.gz").exists()) {
			originalColMap = DoubleMatrixDataset.loadIdentifiers(fileName + ".cols.txt.gz");
		} else {
			throw new FileNotFoundException("File not found: " + fileName + ".cols.txt or " + fileName + ".cols.txt.gz");
		}

	}

	/**
	 * @param rowsToView
	 * @return subset of rows in order of rowsToView
	 * @throws IOException
	 */
	public DoubleMatrixDataset<String, String> loadSubsetOfRowsBinaryDoubleData(LinkedHashSet<String> rowsToView) throws IOException {
		return DoubleMatrixDataset.loadSubsetOfRowsBinaryDoubleData(fileName, rowsToView, originalRowMap, originalColMap);
	}

	/**
	 * @param rowsToView
	 * @return subset of rows in order of rowsToView
	 * @throws IOException
	 */
	public DoubleMatrixDataset<String, String> loadSubsetOfRowsBinaryDoubleData(String[] rowsToView) throws IOException, Exception {
		return loadSubsetOfRowsBinaryDoubleData(Arrays.asList(rowsToView));
	}

	/**
	 * @param rowsToView
	 * @return subset of rows in order of rowsToView
	 * @throws IOException
	 */
	public DoubleMatrixDataset<String, String> loadSubsetOfRowsBinaryDoubleData(Collection<String> rowsToView) throws IOException, Exception {

		LinkedHashSet<String> rowsToViewHash = new LinkedHashSet<>(rowsToView);

		if (rowsToViewHash.size() != rowsToView.size()) {

			StringBuilder duplicateRowsRequested = new StringBuilder();

			HashSet<String> rowsSeen = new HashSet<>();

			for (String rowToView : rowsToView) {

				if (!rowsSeen.add(rowToView)) {
					duplicateRowsRequested.append(rowToView);
					duplicateRowsRequested.append(";");
				}

			}

			throw new Exception("Duplicates in rows not allowed. Requested duplicate values: " + duplicateRowsRequested);
		}

		return loadSubsetOfRowsBinaryDoubleData(rowsToViewHash);

	}

	public int totalRowsInDataset() {
		return originalRowMap.size();
	}

	public int totalColsInDataset() {
		return originalColMap.size();
	}

	public Map<String, Integer> getOriginalRowMap() {
		return Collections.unmodifiableMap(originalRowMap);
	}

	public Map<String, Integer> getOriginalColMap() {
		return Collections.unmodifiableMap(originalColMap);
	}


}
