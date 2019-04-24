/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.matrix2;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Set;
import static umcg.genetica.math.matrix2.DoubleMatrixDataset.*;

/**
 *
 * @author patri
 */
public class DoubleMatrixDatasetRowIterable implements Iterable<double[]> {

	private final int nrRows;
	private final int nrCols;
	private final RandomAccessFile in;
	private final LinkedHashMap<String, Integer> originalRowMap;
	private final LinkedHashMap<String, Integer> originalColMap;

	public DoubleMatrixDatasetRowIterable(String fileName) throws IOException {

		//Now load the row and column identifiers from files
		originalRowMap = loadIdentifiers(fileName + ".rows.txt");
		originalColMap = loadIdentifiers(fileName + ".cols.txt");

		final File fileBinary = new File(fileName + ".dat");
		in = new RandomAccessFile(fileBinary, "r");

		byte[] bytes = new byte[4];
		in.read(bytes, 0, 4);
		nrRows = byteArrayToInt(bytes);
		in.read(bytes, 0, 4);
		nrCols = byteArrayToInt(bytes);

		if (nrRows != originalRowMap.size()) {
			throw new RuntimeException("Matrix at: " + fileName + " does not have expected number of rows");
		}

		if (nrCols != originalColMap.size()) {
			throw new RuntimeException("Matrix at: " + fileName + " does not have expected number of cols");
		}

	}

	public int getNrRows() {
		return nrRows;
	}

	public int getNrCols() {
		return nrCols;
	}

	/**
	 * Backed by linked set so will be ordered
	 *
	 * @return
	 */
	public Set<String> getRows() {
		return Collections.unmodifiableSet(originalRowMap.keySet());
	}

	/**
	 * Backed by linked set so will be ordered
	 *
	 * @return
	 */
	public Set<String> getCols() {
		return Collections.unmodifiableSet(originalColMap.keySet());
	}

	@Override
	public Iterator<double[]> iterator() {
		return new DoubleMatrixDatasetRowIterator();
	}

	private class DoubleMatrixDatasetRowIterator implements Iterator<double[]> {

		private final byte[] buffer;
		private int nextRow;
		private long bits;

		public DoubleMatrixDatasetRowIterator() {
			buffer = new byte[nrCols * 8];
			nextRow = 0;
		}

		@Override
		public boolean hasNext() {
			return nextRow < nrRows;
		}

		@Override
		public double[] next() {

			try {
				final double[] result = new double[nrCols];
				in.read(buffer, 0, nrCols * 8);
				int bufferLoc = 0;
				for (int col = 0; col < nrCols; col++) {
					bits = (long) (buffer[bufferLoc++]) << 56
							| (long) (0xff & buffer[bufferLoc++]) << 48
							| (long) (0xff & buffer[bufferLoc++]) << 40
							| (long) (0xff & buffer[bufferLoc++]) << 32
							| (long) (0xff & buffer[bufferLoc++]) << 24
							| (long) (0xff & buffer[bufferLoc++]) << 16
							| (long) (0xff & buffer[bufferLoc++]) << 8
							| (long) (0xff & buffer[bufferLoc++]);

					result[col] = Double.longBitsToDouble(bits);
				}
				nextRow++;

				return result;

			} catch (IOException ex) {
				throw new RuntimeException(ex);
			}
		}
	}

}
