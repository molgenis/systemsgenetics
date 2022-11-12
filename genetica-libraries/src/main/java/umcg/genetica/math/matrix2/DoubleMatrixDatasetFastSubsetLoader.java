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
import java.util.Set;

/**
 * @author patri
 */
public class DoubleMatrixDatasetFastSubsetLoader {

	private final LinkedHashMap<String, Integer> originalRowMap;
	private final LinkedHashMap<String, Integer> originalColMap;
	private final String fileName;
	private final DoubleMatrixDatasetRowCompressedReader matrixReader;

	public DoubleMatrixDatasetFastSubsetLoader(String fileName) throws IOException {
		this.fileName = fileName;

		if (fileName.endsWith(".datg") || new File(fileName + ".datg").exists()) {
			matrixReader = new DoubleMatrixDatasetRowCompressedReader(fileName);
			originalRowMap = null;
			originalColMap = null;
		} else {
			matrixReader = null;
			// strip .dat or .dat.gz from filename
			if (fileName.endsWith(".dat")) {
				fileName = fileName.substring(0, fileName.length() - 4);
			} else if (fileName.endsWith(".dat.gz")) {
				fileName = fileName.substring(0, fileName.length() - 7);
			}

			//Now load the row and column identifiers from files
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
	}

	/**
	 * @param rowsToView
	 * @return subset of rows in order of rowsToView
	 * @throws IOException
	 */
	public DoubleMatrixDataset<String, String> loadSubsetOfRowsBinaryDoubleData(LinkedHashSet<String> rowsToView) throws IOException {
		if (matrixReader == null) {
			return DoubleMatrixDataset.loadSubsetOfRowsBinaryDoubleData(fileName, rowsToView, originalRowMap, originalColMap);
		} else {
			return matrixReader.loadSubsetOfRows(rowsToView);
		}
	}

	/**
	 * @param rowsToView
	 * @return subset of rows in order of rowsToView
	 * @throws IOException
	 */
	public DoubleMatrixDataset<String, String> loadSubsetOfRowsBinaryDoubleData(String[] rowsToView) throws IOException, Exception {
		if (matrixReader == null) {
			return loadSubsetOfRowsBinaryDoubleData(Arrays.asList(rowsToView));
		} else {
			return matrixReader.loadSubsetOfRows(rowsToView);
		}
	}

	/**
	 * @param rowsToView
	 * @return subset of rows in order of rowsToView
	 * @throws IOException
	 */
	public DoubleMatrixDataset<String, String> loadSubsetOfRowsBinaryDoubleData(Collection<String> rowsToView) throws IOException, Exception {

		if (matrixReader == null) {
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
		} else {
			return matrixReader.loadSubsetOfRows(rowsToView);
		}

	}

	public int totalRowsInDataset() {
		if (matrixReader == null) {
			return originalRowMap.size();
		} else {
			return matrixReader.getNumberOfRows();
		}

	}

	public int totalColsInDataset() {
		if (matrixReader == null) {
			return originalColMap.size();
		} else {
			return matrixReader.getNumberOfColumns();
		}

	}

	public Set<String> getOriginalRowMap() {
		if (matrixReader == null) {
			return Collections.unmodifiableSet(originalRowMap.keySet());
		} else {
			return matrixReader.getRowIdentifiers();
		}

	}

	public Set<String> getOriginalColMap() {
		if (matrixReader == null) {
			return Collections.unmodifiableSet(originalColMap.keySet());
		} else {
			return matrixReader.getColumnIdentifiers();
		}
	}

}
