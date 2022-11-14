/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.matrix2;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Set;

import static umcg.genetica.math.matrix2.DoubleMatrixDataset.*;

/**
 * @author patri
 */
public class DoubleMatrixDatasetRowIterable implements Iterable<double[]> {

	private final int nrRows;
	private final int nrCols;
	private final RandomAccessFile in;
	private final DoubleMatrixDatasetRowCompressedReader matrixReader;
	private final LinkedHashMap<String, Integer> originalRowMap;
	private final LinkedHashMap<String, Integer> originalColMap;

	public DoubleMatrixDatasetRowIterable(String fileName) throws IOException {

		if (fileName.endsWith(".datg") || new File(fileName + ".datg").exists()) {
			
			in = null;//not used by new format
			
			matrixReader = new DoubleMatrixDatasetRowCompressedReader(fileName);
			
			nrRows = matrixReader.getNumberOfRows();
			nrCols = matrixReader.getNumberOfColumns();
			
			originalRowMap = null;
			originalColMap = null;
			
		} else {
			
			//below is support for the old .dat matrix format

			matrixReader = null;//not used by old
			
			File fileBinary = null;
			if (fileName.endsWith(".dat") || fileName.endsWith(".dat.gz")) {
				fileBinary = new File(fileName);
				if (fileName.endsWith(".dat")) {
					fileName = fileName.substring(0, fileName.length() - 4);
				} else if (fileName.endsWith(".dat.gz")) {
					fileName = fileName.substring(0, fileName.length() - 7);
				}
			} else {
				if (new File(fileName + ".dat").exists()) {
					fileBinary = new File(fileName + ".dat");
				} else if (new File(fileName + ".dat.gz").exists()) {
					fileBinary = new File(fileName + ".dat.gz");
				}
			}

			if (fileBinary == null || !fileBinary.exists()) {
				throw new FileNotFoundException("File not found: " + fileName + ".dat or " + fileName + ".dat.gz");
			}

			//Now load the row and column identifiers from files
			if (new File(fileName + ".rows.txt").exists()) {
				originalRowMap = loadIdentifiers(fileName + ".rows.txt");
			} else if (new File(fileName + ".rows.txt.gz").exists()) {
				originalRowMap = loadIdentifiers(fileName + ".rows.txt.gz");
			} else {
				throw new FileNotFoundException("File not found: " + fileName + ".rows.txt or " + fileName + ".rows.txt.gz");
			}

			if (new File(fileName + ".cols.txt").exists()) {
				originalColMap = loadIdentifiers(fileName + ".cols.txt");
			} else if (new File(fileName + ".cols.txt.gz").exists()) {
				originalColMap = loadIdentifiers(fileName + ".cols.txt.gz");
			} else {
				throw new FileNotFoundException("File not found: " + fileName + ".cols.txt or " + fileName + ".cols.txt.gz");
			}

			if (fileBinary.getName().endsWith(".dat.gz")) {
				// move to TMP directory first
				System.out.println("Attempting random access on gzipped matrix: " + fileBinary.getName());
				fileBinary = DoubleMatrixDataset.unzipToTMP(fileBinary);
			}

			in = new RandomAccessFile(fileBinary, "r");

			byte[] bytes = new byte[4];
			in.read(bytes, 0, 4);
			nrRows = byteArrayToInt(bytes);
			in.read(bytes, 0, 4);
			nrCols = byteArrayToInt(bytes);

			if (nrRows != originalRowMap.size()) {
				throw new RuntimeException("Matrix at: " + fileName + " does not have expected number of rows. Expected " + originalRowMap.size() + " but found " + nrRows);
			}

			if (nrCols != originalColMap.size()) {
				throw new RuntimeException("Matrix at: " + fileName + " does not have expected number of cols. Expected " + originalColMap.size() + " but found " + nrCols);
			}

			System.out.println(fileName + " is " + nrRows + " rows and " + nrCols + " cols");
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
		if(matrixReader == null){
			return Collections.unmodifiableSet(originalRowMap.keySet());
		} else {
			return matrixReader.getRowIdentifiers();
		}
		
	}

	/**
	 * Backed by linked set so will be ordered
	 *
	 * @return
	 */
	public Set<String> getCols() {
		if(matrixReader == null){
			return Collections.unmodifiableSet(originalColMap.keySet());
		} else {
			return matrixReader.getColumnIdentifiers();
		}		
	}

	@Override
	public Iterator<double[]> iterator() {
		if(matrixReader == null){
			return new DoubleMatrixDatasetRowIterator(in);
		} else {
			return new DoubleMatrixDatasetRowIteratorDatg(matrixReader);
		}
		
	}

	private class DoubleMatrixDatasetRowIterator implements Iterator<double[]> {

		private final byte[] buffer;
		private int nextRow;
		private long bits;
		private final RandomAccessFile in;

		public DoubleMatrixDatasetRowIterator(RandomAccessFile in) {
			buffer = new byte[nrCols * 8];
			nextRow = 0;
			this.in = in;
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
	
	private class DoubleMatrixDatasetRowIteratorDatg implements Iterator<double[]> {

		private int nextRow;
		private final DoubleMatrixDatasetRowCompressedReader matrixReader;

		public DoubleMatrixDatasetRowIteratorDatg(DoubleMatrixDatasetRowCompressedReader matrixReader) {
			this.matrixReader = matrixReader;
			nextRow = 0;
		}

		@Override
		public boolean hasNext() {
			return nextRow < nrRows;
		}

		@Override
		public double[] next() {

			try {
				return matrixReader.loadSingleRow(nextRow++);
			} catch (IOException ex) {
				throw new RuntimeException(ex);
			}
			
		}
	}

	public void close() throws IOException {
		in.close();
	}

}
