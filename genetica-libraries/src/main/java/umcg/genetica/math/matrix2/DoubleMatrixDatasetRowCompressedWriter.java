/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.matrix2;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;
import com.opencsv.CSVWriter;
import gnu.trove.list.array.TLongArrayList;
import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.time.Instant;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPOutputStream;
import net.jpountz.lz4.LZ4BlockOutputStream;
import org.apache.commons.io.output.CountingOutputStream;

/**
 *
 * @author patri
 */
public class DoubleMatrixDatasetRowCompressedWriter {

	protected static final byte[] MAGIC_BYTES = {85, 77, 67, 71};
	protected static final long APPENDIX_BYTE_LENGTH
			= 8 //timestamp
			+ 4 //number rows
			+ 4 //number columns
			+ 8 //start index block
			+ 8 //start meta data block
			+ 4 //rows per block
			+ 4 //reserve
			+ 4; //magic bytes

	private final CSVWriter rowNamesWriter;
	private final TLongArrayList blockIndices;
	private final String[] outputLine = new String[1];//used for CSV writer
	private final CountingOutputStream matrixFileWriter;//use counting to check how much byte each row used
	private final int numberOfColumns;
	private final int rowsPerBlock;
	private int rowsInCurrentBlock = -1;// -1 indicaties a new block is needed should more rows be added
	private LZ4BlockOutputStream blockCompressionWriter;
	private final byte rowBuffer[];
	private final int bytesPerRow;
	private final String datasetName;
	private final String dataOnRows;//For instance genes
	private final String dataOnCols;//For instance pathays
	private int numberOfRows = 0;
	private final int blockSize;
	private final static int MAX_ROWS = Integer.MAX_VALUE - 8;

	public DoubleMatrixDatasetRowCompressedWriter(String path, final List<Object> columns) throws FileNotFoundException, IOException {
		this(path, columns, "", "", "");
	}

	public DoubleMatrixDatasetRowCompressedWriter(String path, final List<Object> columns, int rowsPerBlock) throws FileNotFoundException, IOException {
		this(path, columns, rowsPerBlock, "", "", "");
	}

	public DoubleMatrixDatasetRowCompressedWriter(String path, final List<Object> columns, String datasetName, String dataOnRows, String dataOnCols) throws FileNotFoundException, IOException {
		this(path, columns, columns.size() < 100 ? (99 + columns.size()) / columns.size() : 1, datasetName, dataOnRows, dataOnCols);
		//(x + y - 1) / y; = ceiling int division. -1 is hard coded in 99
	}

	public DoubleMatrixDatasetRowCompressedWriter(String path, final List<Object> columns, int rowsPerBlock, String datasetName, String dataOnRows, String dataOnCols) throws FileNotFoundException, IOException {

		if ((columns.size() * 8l) > 33554432) {
			//this limit is the max block size of lz4 blocks. Can be solved by using multiple block per row but that is not implemented
			throw new IOException("Too many columns to write " + columns.size() + " max is: " + 33554432 / 8);
		}

		if ((columns.size() * 8l * rowsPerBlock) > 33554432) {
			throw new IOException("Too many rows block");
		}

		this.rowsPerBlock = rowsPerBlock;
		this.bytesPerRow = columns.size() * 8;
		rowBuffer = new byte[this.bytesPerRow];

		this.datasetName = datasetName;
		this.dataOnRows = dataOnRows;
		this.dataOnCols = dataOnCols;

		if (path.endsWith(".datg")) {
			path = path.substring(0, path.length() - 5);
		} else if (path.endsWith(".dat")) {
			path = path.substring(0, path.length() - 4);
		}

		final File matrixFile = new File(path + ".datg");
		final File rowFile = new File(path + ".rows.txt.gz");
		final File colFile = new File(path + ".cols.txt.gz");

		numberOfColumns = columns.size();
		final CSVWriter colNamesWriter = new CSVWriter(new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(colFile)), "UTF8")), '\t', '\0', '\0', "\n");
		for (Object col : columns) {
			outputLine[0] = col.toString();
			colNamesWriter.writeNext(outputLine);
		}
		colNamesWriter.close();

		rowNamesWriter = new CSVWriter(new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(rowFile)), "UTF8")), '\t', '\0', '\0', "\n");

		blockIndices = new TLongArrayList(1000);

		matrixFileWriter = new CountingOutputStream(new BufferedOutputStream(new FileOutputStream(matrixFile), 262144));

		blockSize = bytesPerRow * rowsPerBlock;

	}

	/**
	 * Write double array to a row. Can be made faster with dedicated
	 * implementation instead of converting to DoubleMatrix1D
	 *
	 * @param rowName
	 * @param rowData
	 * @throws IOException
	 */
	public final void addRow(final String rowName, final double[] rowData) throws IOException {
		addRow(rowName, new DenseDoubleMatrix1D(rowData));
	}

	public synchronized final void addRow(final String rowName, final DoubleMatrix1D rowData) throws IOException {

		if (++numberOfRows >= MAX_ROWS) {
			throw new IOException("Reached max number of rows: " + numberOfRows);
		}

		outputLine[0] = rowName;
		rowNamesWriter.writeNext(outputLine);

		if (rowsInCurrentBlock == -1) {
			startBlock();
		}

		int b = 0;
		for (int c = 0; c < numberOfColumns; ++c) {
			long v = Double.doubleToLongBits(rowData.getQuick(c));
			rowBuffer[b++] = (byte) (v >>> 56);
			rowBuffer[b++] = (byte) (v >>> 48);
			rowBuffer[b++] = (byte) (v >>> 40);
			rowBuffer[b++] = (byte) (v >>> 32);
			rowBuffer[b++] = (byte) (v >>> 24);
			rowBuffer[b++] = (byte) (v >>> 16);
			rowBuffer[b++] = (byte) (v >>> 8);
			rowBuffer[b++] = (byte) (v);

		}
		blockCompressionWriter.write(rowBuffer, 0, bytesPerRow);

		if (++rowsInCurrentBlock >= rowsPerBlock) {
			closeBlock();
			//don't start new block because unknown if more data will come
		}

	}

	private void startBlock() {
		blockIndices.add(matrixFileWriter.getByteCount());//Set the start byte for this block
		blockCompressionWriter = new LZ4BlockOutputStream(matrixFileWriter, blockSize);
		rowsInCurrentBlock = 0;
	}

	private void closeBlock() throws IOException {
		blockCompressionWriter.finish();
		rowsInCurrentBlock = -1;
	}

	public synchronized final void close() throws IOException {
		//here write row indicies to end of file
		//Then number of row and columns
		//Finaly start of row indices end block

		if (rowsInCurrentBlock > 0) {//>0 is correct this is count not index
			closeBlock();
		}

		final long startOfIndexBlock = matrixFileWriter.getByteCount();

		final LZ4BlockOutputStream rowIndicesCompression = new LZ4BlockOutputStream(matrixFileWriter);
		final DataOutputStream rowIndicesWriter = new DataOutputStream(rowIndicesCompression);

		final int numberOfBlocks = blockIndices.size();
		for (int r = 0; r < numberOfBlocks; ++r) {
			rowIndicesWriter.writeLong(blockIndices.get(r));
		}
		rowIndicesCompression.finish();

		final DataOutputStream metaDataBlockWriter = new DataOutputStream(matrixFileWriter);

		final long startOfMetaDataBlock = matrixFileWriter.getByteCount();

		metaDataBlockWriter.writeUTF(datasetName);
		metaDataBlockWriter.writeUTF(dataOnRows);
		metaDataBlockWriter.writeUTF(dataOnCols);

		metaDataBlockWriter.writeLong(Instant.now().getEpochSecond());
		metaDataBlockWriter.writeInt(numberOfRows);
		metaDataBlockWriter.writeInt(numberOfColumns);
		metaDataBlockWriter.writeLong(startOfIndexBlock);
		metaDataBlockWriter.writeLong(startOfMetaDataBlock);
		metaDataBlockWriter.writeInt(rowsPerBlock);
		metaDataBlockWriter.writeInt(0);//reserve
		metaDataBlockWriter.write(MAGIC_BYTES);

		metaDataBlockWriter.close();

		rowNamesWriter.close();
	}

	public static final void saveDataset(final String path, final DoubleMatrixDataset dataset) throws FileNotFoundException, IOException {
		saveDataset(path, dataset, "", "", "");
	}

	public static final void saveDataset(final String path, final DoubleMatrixDataset dataset, final String datasetName, final String dataOnRows, final String dataOnCols) throws FileNotFoundException, IOException {

		final DoubleMatrixDatasetRowCompressedWriter writer = new DoubleMatrixDatasetRowCompressedWriter(path, dataset.getColObjects(), datasetName, dataOnRows, dataOnCols);

		final ArrayList rowNames = dataset.getRowObjects();

		final int totalRows = dataset.rows();
		for (int r = 0; r < totalRows; ++r) {
			writer.addRow(rowNames.get(r).toString(), dataset.getRow(r));
		}

		writer.close();

	}

}
