/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.matrix2;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseLargeDoubleMatrix2D;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.NoSuchElementException;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author MarcJan, Juha, Harm-Jan, Patrick
 * @param <R>
 * @param <C>
 */
public class DoubleMatrixDataset<R extends Comparable, C extends Comparable> {

	static final IOException doubleMatrixDatasetNonUniqueHeaderException = new IOException("Tried to use a non-unique header set in an identifier HashMap");
	static final Logger LOGGER = Logger.getLogger(DoubleMatrixDataset.class.getName());

	public static DoubleMatrixDataset<String, String> loadDoubleTextData(String expressionDataPath, char c) {
		throw new UnsupportedOperationException("Not supported yet.");
	}
	protected DoubleMatrix2D matrix;
	protected LinkedHashMap<R, Integer> hashRows;
	protected LinkedHashMap<C, Integer> hashCols;

	public DoubleMatrixDataset() {
		hashRows = new LinkedHashMap<R, Integer>();
		hashCols = new LinkedHashMap<C, Integer>();
	}

	public DoubleMatrixDataset(int rows, int columns) {
		hashRows = new LinkedHashMap<R, Integer>((int) Math.ceil(rows / 0.75));
		hashCols = new LinkedHashMap<C, Integer>((int) Math.ceil(columns / 0.75));
		if ((rows * (long) columns) < (Integer.MAX_VALUE - 2)) {
			matrix = new DenseDoubleMatrix2D(rows, columns);
		} else {
			matrix = new DenseLargeDoubleMatrix2D(rows, columns);
		}
	}

	public DoubleMatrixDataset(LinkedHashMap<R, Integer> hashRows, LinkedHashMap<C, Integer> hashCols) {
		this.hashRows = hashRows;
		this.hashCols = hashCols;
		if ((hashRows.size() * (long) hashCols.size()) < (Integer.MAX_VALUE - 2)) {
			matrix = new DenseDoubleMatrix2D(hashRows.size(), hashCols.size());
		} else {
			matrix = new DenseLargeDoubleMatrix2D(hashRows.size(), hashCols.size());
		}
	}

	public DoubleMatrixDataset(DoubleMatrix2D matrix, LinkedHashMap<R, Integer> hashRows, LinkedHashMap<C, Integer> hashCols) {
		this.hashRows = hashRows;
		this.hashCols = hashCols;
		this.matrix = matrix;
	}

	public DoubleMatrixDataset(Collection<R> rowNames, Collection<C> colNames) {

		hashRows = new LinkedHashMap<R, Integer>(rowNames.size());
		hashCols = new LinkedHashMap<C, Integer>(colNames.size());

		int i = 0;
		for (R row : rowNames) {
			hashRows.put(row, i);
			++i;
		}

		i = 0;
		for (C col : colNames) {
			hashCols.put(col, i);
			++i;
		}

		if ((hashRows.size() * (long) hashCols.size()) < (Integer.MAX_VALUE - 2)) {
			matrix = new DenseDoubleMatrix2D(hashRows.size(), hashCols.size());
		} else {
			matrix = new DenseLargeDoubleMatrix2D(hashRows.size(), hashCols.size());
		}

	}

	public static DoubleMatrixDataset<String, String> loadDoubleData(String fileName) throws IOException {
		if ((fileName.endsWith(".txt") || fileName.endsWith(".tsv") || fileName.endsWith(".txt.gz"))) {
			return loadDoubleTextData(fileName, "\t");
		} else if (fileName.endsWith(".binary")) {
			return loadDoubleBinaryData(fileName);
		} else {
			throw new IllegalArgumentException("File type must be \".txt\", \".tsv\" or \".txt.gz\" when delimiter is set to: \"tab\" \n Input filename: " + fileName);
		}
	}

	public static DoubleMatrixDataset<String, String> loadDoubleTextData(String fileName, String delimiter) throws IOException {
		if (!(fileName.endsWith(".txt") || fileName.endsWith(".tsv") || fileName.endsWith(".txt.gz"))) {
			throw new IllegalArgumentException("File type must be \".txt\", \".tsv\" or \".txt.gz\" when delimiter is set. \n Input filename: " + fileName);
		}

		Pattern splitPatern = Pattern.compile(delimiter);

		int columnOffset = 1;

		TextFile in = new TextFile(fileName, TextFile.R);
		String str = in.readLine(); // header
		String[] data = splitPatern.split(str);

		int tmpCols = (data.length - columnOffset);

		LinkedHashMap<String, Integer> colMap = new LinkedHashMap<String, Integer>((int) Math.ceil(tmpCols / 0.75));

		for (int s = 0; s < tmpCols; s++) {
			String colName = data[s + columnOffset];
			if (!colMap.containsKey(colName)) {
				colMap.put(colName, s);
			} else {
				LOGGER.warning("Duplicated column name!");
				throw (doubleMatrixDatasetNonUniqueHeaderException);
			}
		}

		int tmpRows = 0;

		while (in.readLine() != null) {
			tmpRows++;
		}
		in.close();

		LinkedHashMap<String, Integer> rowMap = new LinkedHashMap<String, Integer>((int) Math.ceil(tmpRows / 0.75));
		DoubleMatrix2D tmpMatrix;

		if ((tmpRows * (long) tmpCols) < (Integer.MAX_VALUE - 2)) {
			tmpMatrix = new DenseDoubleMatrix2D(tmpRows, tmpCols);
		} else {
			tmpMatrix = new DenseLargeDoubleMatrix2D(tmpRows, tmpCols);
		}
		in.open();
		in.readLine(); // read header
		int row = 0;

		boolean correctData = true;
		while ((str = in.readLine()) != null) {
			data = splitPatern.split(str);

			if (!rowMap.containsKey(data[0])) {
				rowMap.put(data[0], row);
				for (int s = 0; s < tmpCols; s++) {
					double d;
					try {
						d = Double.parseDouble(data[s + columnOffset]);
					} catch (NumberFormatException e) {
						correctData = false;
						d = Double.NaN;
					}
					tmpMatrix.setQuick(row, s, d);
				}
				row++;
			} else {
				LOGGER.warning("Duplicated row name!");
				throw (doubleMatrixDatasetNonUniqueHeaderException);
			}

		}
		if (!correctData) {
			LOGGER.warning("Your data contains NaN/unparseable values!");
		}
		in.close();

		DoubleMatrixDataset<String, String> dataset = new DoubleMatrixDataset<String, String>(tmpMatrix, rowMap, colMap);

		LOGGER.log(Level.INFO, "''{0}'' has been loaded, nrRows: {1} nrCols: {2}", new Object[]{fileName, dataset.matrix.rows(), dataset.matrix.columns()});
		return dataset;
	}

	public static DoubleMatrixDataset<String, String> loadSubsetOfTextDoubleData(String fileName, String delimiter, HashSet<String> desiredRows, HashSet<String> desiredCols) throws IOException {
		if (!(fileName.endsWith(".txt") || fileName.endsWith(".txt.gz") || fileName.endsWith(".tsv") || fileName.endsWith(".tsv.gz"))) {
			throw new IllegalArgumentException("File type must be .txt or .tsv when delimiter is given (given filename: " + fileName + ")");
		}

		LinkedHashSet<Integer> desiredColPos = new LinkedHashSet<Integer>();

		Pattern splitPatern = Pattern.compile(delimiter);

		int columnOffset = 1;

		TextFile in = new TextFile(fileName, TextFile.R);
		String str = in.readLine(); // header
		String[] data = splitPatern.split(str);

		int tmpCols = (data.length - columnOffset);

		LinkedHashMap<String, Integer> colMap = new LinkedHashMap<String, Integer>((int) Math.ceil(tmpCols / 0.75));

		int storedCols = 0;
		for (int s = 0; s < tmpCols; s++) {
			String colName = data[s + columnOffset];
			if (!colMap.containsKey(colName) && (desiredCols == null || desiredCols.contains(colName) || desiredCols.isEmpty())) {
				colMap.put(colName, storedCols);
				desiredColPos.add((s));
				storedCols++;
			} else if (colMap.containsKey(colName)) {
				LOGGER.warning("Duplicated column name!");
				System.out.println("Tried to add: " + colName);
				throw (doubleMatrixDatasetNonUniqueHeaderException);
			}
		}

		LinkedHashSet<Integer> desiredRowPos = new LinkedHashSet<Integer>();
		int rowsToStore = 0;
		int totalRows = 0;
		//System.out.println(desiredRows.toString());
		while ((str = in.readLine()) != null) {
			String[] info = splitPatern.split(str);
			if (desiredRows == null || desiredRows.contains(info[0]) || desiredRows.isEmpty()) {
				rowsToStore++;
				desiredRowPos.add(totalRows);
			}
			totalRows++;
		}
		in.close();

		DoubleMatrix2D matrix;
		if ((rowsToStore * (long) tmpCols) < (Integer.MAX_VALUE - 2)) {
			matrix = new DenseDoubleMatrix2D(rowsToStore, storedCols);
		} else {
			matrix = new DenseLargeDoubleMatrix2D(rowsToStore, storedCols);
		}

		in.open();
		in.readLine(); // read header
		int storingRow = 0;
		totalRows = 0;
		LinkedHashMap<String, Integer> rowMap = new LinkedHashMap<String, Integer>((int) Math.ceil(rowsToStore / 0.75));

		boolean correctData = true;
		while ((str = in.readLine()) != null) {

			if (desiredRowPos.contains(totalRows)) {
				data = splitPatern.split(str);
				if (!rowMap.containsKey(data[0])) {
					rowMap.put(data[0], storingRow);
					int storingCol = 0;
					for (int s : desiredColPos) {
						double d;
						try {
							d = Double.parseDouble(data[s + columnOffset]);
						} catch (NumberFormatException e) {
							correctData = false;
							d = Double.NaN;
						}
						matrix.setQuick(storingRow, storingCol, d);
						storingCol++;
					}
					storingRow++;
				} else if (rowMap.containsKey(data[0])) {
					LOGGER.warning("Duplicated row name!");
					System.out.println("Tried to add: " + data[0]);
					throw (doubleMatrixDatasetNonUniqueHeaderException);
				}
			}
			totalRows++;
		}
		if (!correctData) {
			LOGGER.warning("Your data contains NaN/unparseable values!");
		}
		in.close();

		DoubleMatrixDataset<String, String> dataset = new DoubleMatrixDataset<String, String>(matrix, rowMap, colMap);

		LOGGER.log(Level.INFO, "''{0}'' has been loaded, nrRows: {1} nrCols: {2}", new Object[]{fileName, dataset.matrix.rows(), dataset.matrix.columns()});
		return dataset;
	}

	private static DoubleMatrixDataset<String, String> loadDoubleBinaryData(String fileName) throws FileNotFoundException, IOException {
		//First load the raw binary data:
		File fileBinary = new File(fileName + ".dat");
		BufferedInputStream in;
		int nrRows;
		int nrCols;
		in = new BufferedInputStream(new FileInputStream(fileBinary));
		byte[] bytes = new byte[4];
		in.read(bytes, 0, 4);
		nrRows = byteArrayToInt(bytes);
		in.read(bytes, 0, 4);
		nrCols = byteArrayToInt(bytes);

		DoubleMatrix2D matrix;
		if ((nrRows * (long) nrCols) < (Integer.MAX_VALUE - 2)) {
			matrix = new DenseDoubleMatrix2D(nrRows, nrCols);
		} else {
			matrix = new DenseLargeDoubleMatrix2D(nrRows, nrCols);
		}

		//Now load the row and column identifiers from files
		LinkedHashMap<String, Integer> rowMap = loadIdentifiers(fileName + ".rows.txt");
		LinkedHashMap<String, Integer> colMap = loadIdentifiers(fileName + ".cols.txt");

		byte[] buffer = new byte[nrCols * 8];
		long bits;
		for (int row = 0; row < nrRows; row++) {
			in.read(buffer, 0, nrCols * 8);
			int bufferLoc = 0;
			for (int col = 0; col < nrCols; col++) {
				bits = (long) (0xff & buffer[bufferLoc + 7])
						| (long) (0xff & buffer[bufferLoc + 6]) << 8
						| (long) (0xff & buffer[bufferLoc + 5]) << 16
						| (long) (0xff & buffer[bufferLoc + 4]) << 24
						| (long) (0xff & buffer[bufferLoc + 3]) << 32
						| (long) (0xff & buffer[bufferLoc + 2]) << 40
						| (long) (0xff & buffer[bufferLoc + 1]) << 48
						| (long) (buffer[bufferLoc]) << 56;

				matrix.setQuick(row, col, Double.longBitsToDouble(bits));
				bufferLoc += 8;
			}
		}
		in.close();

		DoubleMatrixDataset<String, String> dataset = new DoubleMatrixDataset<String, String>(matrix, rowMap, colMap);
		LOGGER.log(Level.INFO, "Binary file ''{0}'' has been loaded, nrRows: {1} nrCols: {2}", new Object[]{fileName, nrRows, nrCols});

		return dataset;
	}

	private static LinkedHashMap<String, Integer> loadIdentifiers(String filename) throws IOException {
		TextFile tf = new TextFile(filename, false);
		String[] rowsArr = tf.readAsArray();
		tf.close();
		LinkedHashMap<String, Integer> rowMap = new LinkedHashMap<String, Integer>();
		for (String row : rowsArr) {
			rowMap.put(row, rowMap.size());
		}
		return rowMap;
	}

	public void save(File file, String rowDescriptor) throws IOException {
		TextFile out = new TextFile(file, TextFile.W);

		out.append(rowDescriptor);
		for (C col : hashCols.keySet()) {

			out.append('\t');
			out.append(col.toString());
		}
		out.append('\n');
		int r = 0;
		for (R row : hashRows.keySet()) {
			out.append(row.toString());
			for (int c = 0; c < matrix.columns(); c++) {
				out.append('\t');
				out.append(String.valueOf(matrix.getQuick(r, c)));
			}
			out.append('\n');
			++r;
		}
		out.close();
	}

	public void save(String fileName) throws IOException {
		save(new File(fileName), "-");
	}
        public void save(File fileName) throws IOException {
		save(fileName, "-");
	}
        
        public void save(String fileName, String rowDescriptor) throws IOException {
		save(new File(fileName) , rowDescriptor);
	}
        
        public void saveDice(String fileName) throws IOException {
            saveDice(new File(fileName), "-");
        }
        
        public void saveDice(File fileName) throws IOException {
            saveDice(fileName, "-");
        }
        
        public void saveDice(String fileName, String rowDescriptor) throws IOException {
		saveDice(new File(fileName) , rowDescriptor);
	}
        
	public void saveDice(File fileName, String rowDescriptor) throws IOException {
		TextFile out = new TextFile(fileName, TextFile.W);

		out.append(rowDescriptor);
		for (R row : hashRows.keySet()) {
			out.append('\t');
			out.append(row.toString());
		}
		out.append('\n');

		int c = 0;
		for (C col : hashCols.keySet()) {
			out.append(col.toString());
			for (int r = 0; r < matrix.rows(); r++) {

				out.append('\t');
				out.append(String.valueOf(matrix.getQuick(r, c)));
			}
			out.append('\n');
			++c;
		}
		out.close();
	}

	private static byte[] intToByteArray(int value) {
		return new byte[]{(byte) (value >>> 24),
			(byte) (value >>> 16),
			(byte) (value >>> 8),
			(byte) value};
	}

	private static int byteArrayToInt(byte[] b) {
		return (b[0] << 24)
				+ ((b[1] & 0xff) << 16)
				+ ((b[2] & 0xff) << 8)
				+ (b[3] & 0xff);
	}

	//Getters and setters
	public int rows() {
		return matrix.rows();
	}

	public int columns() {
		return matrix.columns();
	}

	public LinkedHashMap<R, Integer> getHashRows() {
		return hashRows;
	}

	public void setHashRows(LinkedHashMap<R, Integer> hashRows) {
		this.hashRows = hashRows;
	}

	public LinkedHashMap<C, Integer> getHashCols() {
		return hashCols;
	}

	public void setHashCols(LinkedHashMap<C, Integer> hashCols) {
		this.hashCols = hashCols;
	}

	public ArrayList<R> getRowObjects() {
		return new ArrayList<R>(hashRows.keySet());
	}

	public void setRowObjects(List<R> arrayList) throws Exception {
		LinkedHashMap<R, Integer> newHashRows = new LinkedHashMap<R, Integer>((int) Math.ceil(arrayList.size() / 0.75));
		int i = 0;
		for (R s : arrayList) {
			if (!newHashRows.containsKey(s)) {
				newHashRows.put(s, i);
			} else {
				System.out.println("Error, new row names contains dupilcates.");
				throw (doubleMatrixDatasetNonUniqueHeaderException);
			}
			i++;
		}

		this.hashRows = newHashRows;
	}

	public ArrayList<C> getColObjects() {
		return new ArrayList<C>(hashCols.keySet());
	}

	public void setColObjects(List<C> arrayList) throws Exception {
		LinkedHashMap<C, Integer> newHashCols = new LinkedHashMap<C, Integer>((int) Math.ceil(arrayList.size() / 0.75));
		int i = 0;
		for (C s : arrayList) {
			if (!newHashCols.containsKey(s)) {
				newHashCols.put(s, i);
			} else {
				System.out.println("Error, new column names contains dupilcates.");
				throw (doubleMatrixDatasetNonUniqueHeaderException);
			}
			i++;
		}
		this.hashCols = newHashCols;
	}

	public DoubleMatrix2D getMatrix() {
		return matrix;
	}

	public void setMatrix(DoubleMatrix2D matrix) {
		this.matrix = matrix;
	}

	public void setMatrix(double[][] matrix) {
		if ((matrix.length * (long) matrix[0].length) < (Integer.MAX_VALUE - 2)) {
			this.matrix = new DenseDoubleMatrix2D(matrix);
		} else {
			this.matrix = new DenseLargeDoubleMatrix2D(matrix.length, matrix[0].length);
			this.matrix.assign(matrix);
		}
	}

	/**
	 * Order columns
	 *
	 */
	public void OrderOnColumnnames() {
		LinkedHashMap<C, Integer> newColHash = new LinkedHashMap<C, Integer>((int) Math.ceil(this.matrix.columns() / 0.75));
		ArrayList<C> names = this.getColObjects();
		Collections.sort(names);

		int pos = 0;
		for (C name : names) {
			newColHash.put(name, pos);
			pos++;
		}
		reorderCols(newColHash);
	}

	/**
	 * Order rows
	 *
	 */
	public void OrderOnRownames() {
		LinkedHashMap<R, Integer> newRowHash = new LinkedHashMap<R, Integer>((int) Math.ceil(this.matrix.rows() / 0.75));
		ArrayList<R> names = this.getRowObjects();
		Collections.sort(names);

		int pos = -1;
		for (R name : names) {
			pos++;
			newRowHash.put(name, pos);
		}
		reorderRows(newRowHash);

	}

	public void reorderRows(LinkedHashMap<R, Integer> mappingIndex) {
		boolean equal = compareHashRows(mappingIndex, this.hashRows);
		if (!equal) {
			DoubleMatrix2D newRawData;
			if ((this.rows() * (long) this.columns()) < (Integer.MAX_VALUE - 2)) {
				newRawData = new DenseDoubleMatrix2D(this.rows(), this.columns());
			} else {
				newRawData = new DenseLargeDoubleMatrix2D(this.rows(), this.columns());
			}

			for (Map.Entry<R, Integer> ent : mappingIndex.entrySet()) {
				int pos = this.getHashRows().get(ent.getKey());
				for (int s = 0; s < this.columns(); ++s) {
					newRawData.set(ent.getValue(), s, this.getMatrix().get(pos, s));
				}
			}
			this.setHashRows(mappingIndex);
			this.setMatrix(newRawData);
		}

	}

	public void reorderCols(LinkedHashMap<C, Integer> mappingIndex) {
		boolean equal = compareHashCols(mappingIndex, this.hashCols);
		if (!equal) {
			DoubleMatrix2D newRawData;
			if ((this.rows() * (long) this.columns()) < (Integer.MAX_VALUE - 2)) {
				newRawData = new DenseDoubleMatrix2D(this.rows(), this.columns());
			} else {
				newRawData = new DenseLargeDoubleMatrix2D(this.rows(), this.columns());
			}

			for (Map.Entry<C, Integer> ent : mappingIndex.entrySet()) {
				int pos = this.getHashCols().get(ent.getKey());
				for (int p = 0; p < this.rows(); ++p) {
					newRawData.set(p, ent.getValue(), this.getMatrix().get(p, pos));
				}
			}

			this.setHashCols(mappingIndex);
			this.setMatrix(newRawData);
		}
	}

	public DoubleMatrixDataset<C, R> viewDice() {
		return new DoubleMatrixDataset<C, R>(matrix.viewDice(), hashCols, hashRows);
	}

	private boolean compareHashCols(LinkedHashMap<C, Integer> mappingIndex, LinkedHashMap<C, Integer> originalHash) {

		for (Entry<C, Integer> entry : mappingIndex.entrySet()) {
			if (entry.getValue() != originalHash.get(entry.getKey())) {
				return false;
			}
		}
		return true;
	}

	private boolean compareHashRows(LinkedHashMap<R, Integer> mappingIndex, LinkedHashMap<R, Integer> originalHash) {

		for (Entry<R, Integer> entry : mappingIndex.entrySet()) {
			if (entry.getValue() != originalHash.get(entry.getKey())) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Set a element of the dataset.
	 *
	 * @param rowName
	 * @param columnName
	 * @param value
	 */
	public void setElement(R rowName, C columnName, double value) {

		Integer row = hashRows.get(rowName);
		Integer column = hashCols.get(columnName);

		if (row != null && column != null) {
			matrix.setQuick(row, column, value);
		} else {
			if (row == null) {
				throw new NoSuchElementException("Row not found: " + rowName.toString());
			} else {
				throw new NoSuchElementException("Column not found: " + columnName.toString());
			}

		}

	}

	/**
	 * Get specific element.
	 *
	 * @param rowName
	 * @param columnName
	 * @return
	 */
	public double getElement(R rowName, C columnName) {

		Integer row = hashRows.get(rowName);
		Integer column = hashCols.get(columnName);

		if (row != null && column != null) {
			return matrix.getQuick(row, column);
		} else {
			if (row == null) {
				throw new NoSuchElementException("Row not found: " + rowName.toString());
			} else {
				throw new NoSuchElementException("Column not found: " + columnName.toString());
			}
		}
	}
	
	public DoubleMatrix1D getRow (R rowName){
		Integer row = hashRows.get(rowName);
		if (row != null){
			return matrix.viewRow(row);
		} else {
			throw new NoSuchElementException("Row not found: " + rowName.toString());
		}
	}

	/**
	 * Get specific element.
	 *
	 * @param row
	 * @param column
	 * @return
	 */
	public double getElement(int row, int column) {

		return matrix.get(row, column);
	}
	
	public boolean containsRow(R rowId){
		return hashRows.containsKey(rowId);
	}
	
	public boolean containsCol(C colId){
		return hashCols.containsKey(colId);
	}
}
