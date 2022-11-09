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
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.zip.Checksum;
import java.util.zip.GZIPInputStream;
import net.jpountz.lz4.LZ4BlockInputStream;
import net.jpountz.lz4.LZ4Factory;
import net.jpountz.lz4.LZ4FastDecompressor;
import net.jpountz.xxhash.XXHashFactory;
import org.apache.commons.io.input.RandomAccessFileInputStream;


/**
 *
 * @author patri
 */
public class DoubleMatrixDatasetRowCompressedReader {

	private final LinkedHashMap<String, Integer> rowMap;
	private final LinkedHashMap<String, Integer> colMap;
	private final int numberOfRows;
	private final int numberOfColumns;
	private final int bytesPerRow;
	private final RandomAccessFile matrixFileReader;
	private final long[] blockIndices;
	private final File matrixFile;
	private final long startOfEndBlock;
	private final int blockSize;
	private final int rowsPerBlock;

	public DoubleMatrixDatasetRowCompressedReader(String path) throws IOException {

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
		
		bytesPerRow = numberOfColumns * 8;
		
		matrixFileReader.seek(startOfEndBlock);
		DataInputStream blockIndicesReader = new DataInputStream(new LZ4BlockInputStream(new RandomAccessFileInputStream(matrixFileReader, false)));
		blockIndices = new long[numberOfRows];
		for (int r = 0; r < numberOfRows; ++r) {
			blockIndices[r] = blockIndicesReader.readLong();
		}
		blockIndicesReader.close();
		
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

	public synchronized DoubleMatrixDataset<String, String> loadSubsetOfRows(LinkedHashSet<String> rowsToView) throws IOException {


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

			matrixFileReader.seek(blockIndices[rowMap.get(rowName)]);

			DataInputStream reader = new DataInputStream(new GZIPInputStream(new BufferedInputStream(new RandomAccessFileInputStream(matrixFileReader, false), 262144)));

			for (int c = 0; c < numberOfColumns; ++c) {
				matrix.setQuick(rowIndex, c, reader.readDouble());
			}

			reader.available();

		}

		return dataset;

	}

	public synchronized DoubleMatrixDataset<String, String> loadFullDataset() throws IOException {

		final RandomAccessFileInputStream inputStream = new RandomAccessFileInputStream(matrixFileReader, false);//new CountingInputStream(new BufferedInputStream(new RandomAccessFileInputStream(matrixFileReader, false),262144));
		
		
		//make clones to make sure the oringal remains intact
		final DoubleMatrixDataset dataset = new DoubleMatrixDataset((LinkedHashMap) rowMap.clone(), (LinkedHashMap) colMap.clone());
		final DoubleMatrix2D matrix = dataset.getMatrix();
		
		final LZ4FastDecompressor decompressor = LZ4Factory.fastestInstance().fastDecompressor();
		final Checksum hasher = XXHashFactory.fastestInstance().newStreamingHash32(0x9747b28c).asChecksum(); //0x9747b28c is the default seed used by default constructor
		
		final byte[] rowData = new byte[numberOfColumns * 8];
		
		int b;
		LZ4BlockInputStream reader;
		for (int r = 0; r < numberOfRows; ++r) {
		
			matrixFileReader.seek(blockIndices[r]);
		
			reader = new LZ4BlockInputStream(inputStream,decompressor, hasher);

			//Benchmarking revealed that DataInputStream is slow when you can manage the exact amount of data needed
			reader.read(rowData, 0, bytesPerRow);

			b = 0;
			for (int c = 0; c < numberOfColumns; ++c) {

				//More effcient than DataInputStream because we do not need to check the buffer.
				matrix.setQuick(r, c,
						Double.longBitsToDouble(((long) rowData[b++] << 56)
								+ ((long) (rowData[b++] & 255) << 48)
								+ ((long) (rowData[b++] & 255) << 40)
								+ ((long) (rowData[b++] & 255) << 32)
								+ ((long) (rowData[b++] & 255) << 24)
								+ ((rowData[b++] & 255) << 16)
								+ ((rowData[b++] & 255) << 8)
								+ ((rowData[b++] & 255)))
				);
			}
			
			
		}

		return dataset;

	}

}
