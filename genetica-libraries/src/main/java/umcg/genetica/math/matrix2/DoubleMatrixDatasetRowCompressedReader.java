/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.matrix2;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.RandomAccessFile;
import java.time.Instant;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;
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
	private final RandomAccessFileInputStream inputStream;
	private final long[] blockIndices;
	private final File matrixFile;
	private final long startOfEndBlock;
	private final int blockSize;
	private final int numberOfBlocks;
	private final int rowsPerBlock;
	private final long timestampCreated;
	private final String datasetName;
	private final String dataOnRows;//For instnace genes
	private final String dataOnCols;//For instance pathays
	private final byte[] blockData;
	private int currentBlock = -1;
	private final LZ4FastDecompressor decompressor = LZ4Factory.fastestInstance().fastDecompressor();
	private final Checksum hasher = XXHashFactory.fastestInstance().newStreamingHash32(0x9747b28c).asChecksum(); //0x9747b28c is the default seed used by default constructor

	public DoubleMatrixDatasetRowCompressedReader(String path) throws IOException {

		if (path.endsWith(".datg")) {
			path = path.substring(0, path.length() - 5);
		}

		matrixFile = new File(path + ".datg");
		final File rowFile = new File(path + ".rows.txt.gz");
		final File colFile = new File(path + ".cols.txt.gz");

		if (!matrixFile.exists()) {
			throw new IOException("Matrix file not found at: " + matrixFile.getAbsolutePath());
		}

		if (!rowFile.exists()) {
			throw new IOException("File with row identifiers not found at: " + rowFile.getAbsolutePath());
		}

		if (!colFile.exists()) {
			throw new IOException("File with col identifiers not found at: " + colFile.getAbsolutePath());
		}

		rowMap = loadIdentifiers(rowFile);
		colMap = loadIdentifiers(colFile);

		matrixFileReader = new RandomAccessFile(matrixFile, "r");
		inputStream = new RandomAccessFileInputStream(matrixFileReader, false);//this is later used by block reader

		matrixFileReader.seek(matrixFileReader.length() - DoubleMatrixDatasetRowCompressedWriter.APPENDIX_BYTE_LENGTH);

		timestampCreated = matrixFileReader.readLong();

		System.out.println("date: " + Instant.ofEpochSecond(timestampCreated).toString());

		numberOfRows = matrixFileReader.readInt();
		numberOfColumns = matrixFileReader.readInt();

		if (numberOfRows != rowMap.size()) {
			throw new IOException("Number of rows in " + rowFile.getAbsolutePath() + " (" + rowMap.size() + ") does not match number of rows found in " + matrixFile.getAbsolutePath() + "(" + numberOfRows + ")");
		}

		if (numberOfColumns != colMap.size()) {
			throw new IOException("Number of columns in " + colFile.getAbsolutePath() + " (" + colMap.size() + ") does not match number of columns found in " + matrixFile.getAbsolutePath() + "(" + numberOfColumns + ")");
		}

		System.out.println("rows " + numberOfRows);
		System.out.println("cols " + numberOfColumns);

		startOfEndBlock = matrixFileReader.readLong();
		final long startOfMetaDataBlock = matrixFileReader.readLong();

		this.rowsPerBlock = matrixFileReader.readInt();

		if (matrixFileReader.readInt() != 0) {
			throw new RuntimeException("Incompatible version of datg file");
		}

		if (matrixFileReader.readByte() != DoubleMatrixDatasetRowCompressedWriter.MAGIC_BYTES[0]
				|| matrixFileReader.readByte() != DoubleMatrixDatasetRowCompressedWriter.MAGIC_BYTES[1]
				|| matrixFileReader.readByte() != DoubleMatrixDatasetRowCompressedWriter.MAGIC_BYTES[2]
				|| matrixFileReader.readByte() != DoubleMatrixDatasetRowCompressedWriter.MAGIC_BYTES[3]) {
			throw new IOException("Error parsing datg file.");
		}

		matrixFileReader.seek(startOfMetaDataBlock);
		this.datasetName = matrixFileReader.readUTF();
		this.dataOnRows = matrixFileReader.readUTF();
		this.dataOnCols = matrixFileReader.readUTF();

		System.out.println("dataset meta: " + this.datasetName + " " + this.dataOnRows + " " + this.dataOnCols);

		bytesPerRow = numberOfColumns * 8;
		blockSize = bytesPerRow * rowsPerBlock;
		blockData = new byte[blockSize];

		System.out.println("rowsPerBlock: " + rowsPerBlock);

		numberOfBlocks = (numberOfRows + rowsPerBlock - 1) / rowsPerBlock; // (x + y - 1) / y; = ceiling divide

		System.out.println("numberOfBlocks " + numberOfBlocks);

		matrixFileReader.seek(startOfEndBlock);
		final DataInputStream blockIndicesReader = new DataInputStream(new LZ4BlockInputStream(new RandomAccessFileInputStream(matrixFileReader, false)));
		blockIndices = new long[numberOfBlocks];
		for (int block = 0; block < numberOfBlocks; ++block) {
			blockIndices[block] = blockIndicesReader.readLong();
		}
		blockIndicesReader.close();

	}

	public synchronized DoubleMatrixDataset<String, String> loadSubsetOfRows(final String[] rowsToView) throws Exception {
		return loadSubsetOfRows(Arrays.asList(rowsToView));
	}

	public synchronized DoubleMatrixDataset<String, String> loadSubsetOfRows(final Collection<String> rowsToView) throws IOException, Exception {

		final LinkedHashSet<String> rowsToViewHash = new LinkedHashSet<>(rowsToView.size());

		for (String rowToView : rowsToView) {
			rowsToViewHash.add(rowToView);
		}

		if (rowsToViewHash.size() != rowsToView.size()) {

			final StringBuilder duplicateRowsRequested = new StringBuilder();

			final HashSet<String> rowsSeen = new HashSet<>();

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

	public synchronized double[] loadSingleRow(final int r) throws IOException {

		if (r >= numberOfRows) {
			throw new IOException("Requested row " + r + "/" + numberOfRows + " in " + matrixFile.getAbsolutePath());
		}

		int b = prepareForRowRead(r);

		final double[] rowContent = new double[numberOfColumns];

		for (int c = 0; c < numberOfColumns; ++c) {

			//More effcient than DataInputStream because we do not need to check the buffer.
			rowContent[c] = Double.longBitsToDouble(((long) blockData[b++] << 56)
					+ ((long) (blockData[b++] & 255) << 48)
					+ ((long) (blockData[b++] & 255) << 40)
					+ ((long) (blockData[b++] & 255) << 32)
					+ ((long) (blockData[b++] & 255) << 24)
					+ ((blockData[b++] & 255) << 16)
					+ ((blockData[b++] & 255) << 8)
					+ ((blockData[b++] & 255)));
		}

		return rowContent;

	}

	public synchronized DoubleMatrixDataset<String, String> loadSubsetOfRows(final LinkedHashSet<String> rowsToView) throws IOException {

		final LinkedHashMap<String, Integer> rowMapSelection = new LinkedHashMap<>();
		for (String rowToView : rowsToView) {
			if (!rowMap.containsKey(rowToView)) {
				throw new RuntimeException("Matrix at: " + matrixFile.getAbsolutePath() + " does not contain this row: " + rowToView);
			}
			rowMapSelection.put(rowToView, rowMapSelection.size());
		}

		final DoubleMatrixDataset dataset = new DoubleMatrixDataset(rowMapSelection, (LinkedHashMap) colMap.clone());
		final DoubleMatrix2D matrix = dataset.getMatrix();

		int b;
		for (Map.Entry<String, Integer> rowToViewEntry : rowMapSelection.entrySet()) {
			final String rowName = rowToViewEntry.getKey();
			final int rowIndex = rowToViewEntry.getValue();

			b = prepareForRowRead(rowMap.get(rowName));

			for (int c = 0; c < numberOfColumns; ++c) {

				//More effcient than DataInputStream because we do not need to check the buffer.
				matrix.setQuick(rowIndex, c,
						Double.longBitsToDouble(((long) blockData[b++] << 56)
								+ ((long) (blockData[b++] & 255) << 48)
								+ ((long) (blockData[b++] & 255) << 40)
								+ ((long) (blockData[b++] & 255) << 32)
								+ ((long) (blockData[b++] & 255) << 24)
								+ ((blockData[b++] & 255) << 16)
								+ ((blockData[b++] & 255) << 8)
								+ ((blockData[b++] & 255)))
				);
			}

		}

		return dataset;

	}

	public final synchronized DoubleMatrixDataset<String, String> loadFullDataset() throws IOException {

		//make clones to make sure the oringal remains intact
		final DoubleMatrixDataset dataset = new DoubleMatrixDataset((LinkedHashMap) rowMap.clone(), (LinkedHashMap) colMap.clone());
		final DoubleMatrix2D matrix = dataset.getMatrix();

		int b;

		for (int r = 0; r < numberOfRows; ++r) {

			b = prepareForRowRead(r);

			for (int c = 0; c < numberOfColumns; ++c) {

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
			}

		}

		return dataset;

	}

	private int prepareForRowRead(final int r) throws IOException {

		final int neededBlock = r / rowsPerBlock;
		if (currentBlock != neededBlock) {
			loadBlockData(neededBlock);
		}
		return ((r % rowsPerBlock) * bytesPerRow);//start of data in current block

	}

	private void loadBlockData(final int block) throws IOException {
		matrixFileReader.seek(blockIndices[block]);
		final LZ4BlockInputStream reader = new LZ4BlockInputStream(inputStream, decompressor, hasher);

		reader.read(blockData, 0, blockSize);

		currentBlock = block;
	}

	public int getNumberOfRows() {
		return numberOfRows;
	}

	public int getNumberOfColumns() {
		return numberOfColumns;
	}

	public File getMatrixFile() {
		return matrixFile;
	}

	public Set<String> getRowIdentifiers() {
		return Collections.unmodifiableSet(rowMap.keySet());
	}

	public Set<String> getColumnIdentifiers() {
		return Collections.unmodifiableSet(colMap.keySet());
	}

	private LinkedHashMap<String, Integer> loadIdentifiers(final File file) throws IOException {

		final CSVParser identifierFileParser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();

		CSVReader reader = new CSVReaderBuilder(new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file)), "UTF8"))).withCSVParser(identifierFileParser).build();

		LinkedHashMap<String, Integer> map = new LinkedHashMap<>();
		String[] inputLine;
		while ((inputLine = reader.readNext()) != null) {
			map.put(inputLine[0], map.size());
		}

		return (map);

	}
}
