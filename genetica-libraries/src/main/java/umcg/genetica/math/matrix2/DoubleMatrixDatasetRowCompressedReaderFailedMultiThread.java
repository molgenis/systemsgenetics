/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.matrix2;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.RandomAccessFile;
import java.lang.Thread.State;
import java.sql.Time;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.GZIPInputStream;
import org.apache.commons.io.input.RandomAccessFileInputStream;
import umontreal.iro.lecuyer.util.TimeUnit;

/**
 *
 * @author patri
 */
public class DoubleMatrixDatasetRowCompressedReaderFailedMultiThread {

	private final LinkedHashMap<String, Integer> rowMap;
	private final LinkedHashMap<String, Integer> colMap;
	private final int numberOfRows;
	private final int numberOfColumns;
	private final RandomAccessFile matrixFileReader;
	private long[] rowIndices = null;
	private final File matrixFile;
	private final long startOfEndBlock;

	public DoubleMatrixDatasetRowCompressedReaderFailedMultiThread(String path) throws IOException {

		if (path.endsWith(".datg")) {
			path = path.substring(0, path.length() - 5);
		}

		matrixFile = new File(path + ".datg");
		final File rowFile = new File(path + ".rows.txt.gz");
		final File colFile = new File(path + ".cols.txt.gz");

		rowMap = DoubleMatrixDataset.loadIdentifiers(rowFile);
		colMap = DoubleMatrixDataset.loadIdentifiers(colFile);

		matrixFileReader = new RandomAccessFile(matrixFile, "r");

		matrixFileReader.seek(matrixFileReader.length() - DoubleMatrixDatasetRowCompressedWriter.APPENDIX_BYTE_LENGTH);

		numberOfRows = matrixFileReader.readInt();
		numberOfColumns = matrixFileReader.readInt();

		System.out.println("rows " + numberOfRows);
		System.out.println("cols " + numberOfColumns);

		startOfEndBlock = matrixFileReader.readLong();

		if (matrixFileReader.readByte() != DoubleMatrixDatasetRowCompressedWriter.MAGIC_BYTES[0]
				|| matrixFileReader.readByte() != DoubleMatrixDatasetRowCompressedWriter.MAGIC_BYTES[1]
				|| matrixFileReader.readByte() != DoubleMatrixDatasetRowCompressedWriter.MAGIC_BYTES[2]
				|| matrixFileReader.readByte() != DoubleMatrixDatasetRowCompressedWriter.MAGIC_BYTES[3]) {
			throw new IOException("Error parsing datg file.");
		}
	}

	public synchronized DoubleMatrixDataset<String, String> loadSubsetOfRows(String[] rowsToView) throws Exception {
		return loadSubsetOfRows(Arrays.asList(rowsToView));
	}

	public synchronized DoubleMatrixDataset<String, String> loadSubsetOfRows(Collection<String> rowsToView) throws IOException, Exception {

		LinkedHashSet<String> rowsToViewHash = new LinkedHashSet<>(rowsToView.size());

		for (String rowToView : rowsToView) {
			rowsToViewHash.add(rowToView);
		}

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

		return loadSubsetOfRows(rowsToViewHash);

	}

	public synchronized void loadRowIndices() throws IOException {

		matrixFileReader.seek(startOfEndBlock);

		DataInputStream endBlockReader = new DataInputStream(new GZIPInputStream(new BufferedInputStream(new RandomAccessFileInputStream(matrixFileReader, false), 262144)));

		rowIndices = new long[numberOfRows];

		for (int r = 0; r < numberOfRows; ++r) {
			rowIndices[r] = endBlockReader.readLong();
		}

		endBlockReader.available();
		endBlockReader.close();

	}

	public synchronized DoubleMatrixDataset<String, String> loadSubsetOfRows(LinkedHashSet<String> rowsToView) throws IOException {

		if (rowIndices == null) {
			//row indices not yet loaded, doing that now.
			loadRowIndices();
		}

		final LinkedHashMap<String, Integer> rowMapSelection = new LinkedHashMap<>();
		for (String rowToView : rowsToView) {
			if (!rowMap.containsKey(rowToView)) {
				throw new RuntimeException("Matrix at: " + matrixFile.getAbsolutePath() + " does not contain this row: " + rowToView);
			}
			rowMapSelection.put(rowToView, rowMapSelection.size());
		}

		DoubleMatrixDataset dataset = new DoubleMatrixDataset(rowMapSelection, (LinkedHashMap) colMap.clone());
		DoubleMatrix2D matrix = dataset.getMatrix();

		for (Map.Entry<String, Integer> rowToViewEntry : rowMapSelection.entrySet()) {
			final String rowName = rowToViewEntry.getKey();
			final int rowIndex = rowToViewEntry.getValue();

			matrixFileReader.seek(rowIndices[rowMap.get(rowName)]);

			DataInputStream reader = new DataInputStream(new GZIPInputStream(new BufferedInputStream(new RandomAccessFileInputStream(matrixFileReader, false), 262144)));

			for (int c = 0; c < numberOfColumns; ++c) {
				matrix.setQuick(rowIndex, c, reader.readDouble());
			}

			reader.available();

		}

		return dataset;

	}

	public synchronized DoubleMatrixDataset<String, String> loadFullDataset() throws IOException, InterruptedException {

		//Move to start of file because we need to read everything
		matrixFileReader.seek(0);

		//Good expercience setting gzip buffer to 2 * number of columns 
		final GZIPInputStream reader = new GZIPInputStream(new BufferedInputStream(new RandomAccessFileInputStream(matrixFileReader, false), 262144));

		//make clones to make sure the oringal remains intact
		final DoubleMatrixDataset dataset = new DoubleMatrixDataset((LinkedHashMap) rowMap.clone(), (LinkedHashMap) colMap.clone());
		final DoubleMatrix2D matrix = dataset.getMatrix();

		//this array is filled by the blockReader thread.
		final byte[] blockData = new byte[numberOfColumns * 8];
		final int bytesPerRow = numberOfColumns * 8;

		BlockReader blockReader = new BlockReader(blockData, bytesPerRow, reader);
		Thread blockReaderThread = new Thread(blockReader);
		blockReaderThread.start();

		int b;
		for (int r = 0; r < numberOfRows; ++r) {

//			if (r % 1000 == 0) {
//				System.out.println("r: " + r);
//			}

//			int j = 0;
//			while (!blockReader.isDone()) {
//				//ensure block reader is not working on something.
//				if (j % 100000 == 0) {
//					System.out.println("stuck waiting for done");
//				}
//			}
			//int j = 0;
//			System.out.println("test0");
//			synchronized (reader) {
//				if (blockReaderThread.getState() != State.WAITING) {
//					System.out.println("test00");
//					reader.wait();
//				}
//			}
//			System.out.println("test1 " + System.currentTimeMillis());
//			synchronized (blockData) {
//				System.out.println("test2 " + System.currentTimeMillis() + blockReaderThread.getState());
//				blockReader.reset();
//				blockData.notify();
//			}
			
			synchronized(blockData){
				blockReader.reset();
				blockData.notify();
			}

			//System.out.println("test r: " + r );
			b = 0;
			for (int c = 0; c < numberOfColumns; ++c) {

				//System.out.println("test");
				int i = 0;
				// >= because of difference in index vs byte count

				synchronized (reader) {
					//System.out.println("test4");
					if ((b + 8) > blockReader.bytesRead) {
						//long x = System.currentTimeMillis();
						reader.wait();
						//System.out.println(System.currentTimeMillis() - x);
					}
				}

				//System.out.println("test r: " + r + " c: " + c + " b: " + b + " read: " + blockReader.bytesRead );
				//More effcient than DataInputStream because we do not need to check the buffer.
				matrix.setQuick(r, c,
						Double.longBitsToDouble(((long) blockData[b++] << 56)
								+ ((long) (blockData[b++] & 255) << 48)
								+ ((long) (blockData[b++] & 255) << 40)
								+ ((long) (blockData[b++] & 255) << 32)
								+ ((long) (blockData[b++] & 255) << 24)
								+ ((blockData[b++] & 255) << 16)
								+ ((blockData[b++] & 255) << 8)
								+ ((blockData[b++] & 255)))
				);

				//System.out.println("test7 " + System.currentTimeMillis());
				synchronized (reader) {
					if (!blockReader.done) {
						reader.wait();
					}
				}
				//System.out.println("test8 " + System.currentTimeMillis());

			}

		}

		reader.available();//this will trigger reading the end of the last gzip block which will check of the checksum of the read data is okay

		reader.close();

		return dataset;

	}

	private class BlockReader implements Runnable {

		private final byte[] blockData;
		private int bytesRead;
		private final int bytesPerBlock;
		private final InputStream reader;
		private boolean done = false;

		public BlockReader(byte[] blockData, int bytesPerRow, InputStream reader) {
			this.blockData = blockData;
			this.bytesPerBlock = bytesPerRow;
			this.reader = reader;
		}

		@Override
		public void run() {

			while (true) {
				try {

					//System.out.println("TestB " + System.currentTimeMillis());
					//System.out.println("processing block");
					//Benchmarking revealed that DataInputStream is slow when you can manage the exact amount of data needed
					//However GZIPInputStream does not return large chunks when reading
					///Therefore small loop to read all data in the row. 
					//This is over 50% faster than wrapping GZIPInputStream in DataInputStream
					int i = 0;
					while (bytesRead < bytesPerBlock) {
						//i++;
						bytesRead += reader.read(blockData, bytesRead, bytesPerBlock - bytesRead);
//						if (i % 100000 == 0) {
//							System.out.println("long processing " + i);
//						}
						//System.out.println("TestC " + System.currentTimeMillis());
						synchronized (reader) {
							reader.notify();
						}
					}

					synchronized (reader) {
						done = true;
						reader.notify();
					}

					//System.out.println("TestC");
					try {
						synchronized (blockData) {
							if(bytesRead ==0 ){
								continue;
							}
							blockData.wait();
						}
					} catch (InterruptedException ex) {
						return;
					}

					//System.out.println("TestD");

					//System.out.println("TestD " + System.currentTimeMillis());
					//System.out.println("done");
				} catch (IOException ex) {
					throw new RuntimeException(ex);
				}
			}

		}

		public void reset() {
			bytesRead = 0;
		}

		public int getBytesRead() {
			return bytesRead;
		}

	}

}
