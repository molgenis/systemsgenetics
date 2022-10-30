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
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPOutputStream;
import org.apache.commons.io.output.CountingOutputStream;

/**
 *
 * @author patri
 */
public class DoubleMatrixDatasetRowCompressedWriter {

	protected static final byte[] MAGIC_BYTES = {85, 77, 67, 71};

	private final CSVWriter rowNamesWriter;
	private final TLongArrayList rowIndices;
	private final String[] outputLine = new String[1];//used for CSV writer
	private final CountingOutputStream matrixFileWriter;//use counting to check how much byte each row used
	private final int numberOfColumns;

	public DoubleMatrixDatasetRowCompressedWriter(String path, List<Object> columns) throws FileNotFoundException, IOException {

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

		rowIndices = new TLongArrayList(1000);

		matrixFileWriter = new CountingOutputStream(new BufferedOutputStream(new FileOutputStream(matrixFile), 262144));
		matrixFileWriter.write(MAGIC_BYTES);

	}

	/**
	 * Write double array to a row. Can be made faster with dedicated
	 * implementation instead of converting to DoubleMatrix1D
	 *
	 * @param rowName
	 * @param rowData
	 * @throws IOException
	 */
	public void addRow(String rowName, double[] rowData) throws IOException {
		addRow(rowName, new DenseDoubleMatrix1D(rowData));
	}

	public void addRow(String rowName, DoubleMatrix1D rowData) throws IOException {

		outputLine[0] = rowName;
		rowNamesWriter.writeNext(outputLine);

		rowIndices.add(matrixFileWriter.getByteCount());//Set the start byte for this row

		final GZIPOutputStream rowCompression = new GZIPOutputStream(matrixFileWriter);
		final DataOutputStream rowWriter = new DataOutputStream(rowCompression);

		for (int c = 0; c < numberOfColumns; ++c) {
			rowWriter.writeDouble(rowData.getQuick(c));
		}

		rowCompression.finish();

	}

	public void close() throws IOException {
		//here write row indicies to end of file
		//Then number of row and columns
		//Finaly start of row indices end block

		final long startOfEndBlock = matrixFileWriter.getByteCount();

		final int numberOfRows = rowIndices.size();
		
		final DataOutputStream endBlockWriter = new DataOutputStream(matrixFileWriter);

		for (int r = 0; r < numberOfRows; ++r) {
			endBlockWriter.writeLong(rowIndices.get(r));
		}
		endBlockWriter.writeInt(numberOfRows);
		endBlockWriter.writeInt(numberOfColumns);
		endBlockWriter.writeLong(startOfEndBlock);

		endBlockWriter.close();
		
		rowNamesWriter.close();
	}

	public static final void saveDataset(String path, DoubleMatrixDataset dataset) throws FileNotFoundException, IOException {

		final DoubleMatrixDatasetRowCompressedWriter writer = new DoubleMatrixDatasetRowCompressedWriter(path, dataset.getColObjects());

		final ArrayList rowNames = dataset.getRowObjects();

		final int totalRows = dataset.rows();
		for (int r = 0; r < totalRows; ++r) {
			writer.addRow(rowNames.get(r).toString(), dataset.getRow(r));
		}

		writer.close();

	}

}
