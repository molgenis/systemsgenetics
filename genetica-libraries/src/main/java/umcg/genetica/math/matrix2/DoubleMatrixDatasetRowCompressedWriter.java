/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.matrix2;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
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
import java.nio.charset.StandardCharsets;
import java.security.DigestOutputStream;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.time.Instant;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.SequencedCollection;
import java.util.zip.GZIPOutputStream;


import org.apache.commons.compress.compressors.lz4.FramedLZ4CompressorOutputStream;
import org.apache.commons.io.output.CountingOutputStream;

/**
 * @author patri
 */
@SuppressWarnings("unused")
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
    private static final int MAX_BLOCK_SIZE = 2_000_000_000;
    private static final int MAX_COLUMNS = 250_000_000;

    private final CSVWriter rowNamesWriter;
    private final TLongArrayList blockIndices;
    private final String[] outputLine = new String[1];//used for CSV writer
    private final CountingOutputStream matrixFileWriter;//use counting to check how much byte each row used
    private final int numberOfColumns;
    private final int rowsPerBlock;
    private int rowsInCurrentBlock = -1;// -1 indicates a new block is needed should more rows be added
    private FramedLZ4CompressorOutputStream blockCompressionWriter;
    private final byte[] rowBuffer;
    private final int bytesPerRow;
    private final String datasetName;
    private final String dataOnRows;//For instance genes
    private final String dataOnCols;//For instance pathways
    private int numberOfRows = 0;
    private final static int MAX_ROWS = Integer.MAX_VALUE - 8;
    private final MessageDigest matrixFileDigest;
    private final MessageDigest rowFileDigest;
    private final MessageDigest colFileDigest;
    private final File matrixFile;
    private final File rowFile;
    private final File colFile;
    private final File checksumFile;
    private final HashSet<String> rowNames = new HashSet<>();

    public DoubleMatrixDatasetRowCompressedWriter(String path, final SequencedCollection<String> columns) throws FileNotFoundException, IOException {
        this(path, columns, "", "", "");
    }

    public DoubleMatrixDatasetRowCompressedWriter(String path, final SequencedCollection<String> columns, int rowsPerBlock) throws FileNotFoundException, IOException {
        this(path, columns, rowsPerBlock, "", "", "");
    }

    public DoubleMatrixDatasetRowCompressedWriter(String path, final SequencedCollection<String> columns, String datasetName, String dataOnRows, String dataOnCols) throws FileNotFoundException, IOException {
        this(path, columns, columns.size() < 100 ? (99 + columns.size()) / columns.size() : 1, datasetName, dataOnRows, dataOnCols);
        //(x + y - 1) / y; = ceiling int division. -1 is hard coded in 99
    }

    public DoubleMatrixDatasetRowCompressedWriter(String path, final SequencedCollection<String> columns, int rowsPerBlock, String datasetName, String dataOnRows, String dataOnCols) throws FileNotFoundException, IOException {

        if (columns.size() > MAX_COLUMNS) {
            //this limit is the max block size of lz4 blocks. Can be solved by using multiple block per row but that is not implemented
            throw new IOException("Too many columns to write " + columns.size() + " max is: " + MAX_COLUMNS);
        }

        if ((columns.size() * 8L * rowsPerBlock) > MAX_BLOCK_SIZE) {
            throw new IOException("Too many rows per block");
        }

        try {
            matrixFileDigest = MessageDigest.getInstance("SHA-256");
            rowFileDigest = MessageDigest.getInstance("SHA-256");
            colFileDigest = MessageDigest.getInstance("SHA-256");
        } catch (NoSuchAlgorithmException e) {
            throw new RuntimeException(e);
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

        matrixFile = new File(path + ".datg");
        rowFile = new File(path + ".rows.txt.gz");
        colFile = new File(path + ".cols.txt.gz");
        checksumFile = new File(path + ".sha256");

        if (matrixFile.getParentFile() != null && !matrixFile.getParentFile().exists()) {
            if (!matrixFile.getParentFile().mkdirs()) {
                throw new IOException("Unable to create directory " + matrixFile.getParent());
            }
        }

        final HashSet<String> colNames = new HashSet<>();
        numberOfColumns = columns.size();
        final CSVWriter colNamesWriter = new CSVWriter(new OutputStreamWriter(new GZIPOutputStream(new DigestOutputStream(new FileOutputStream(colFile), colFileDigest)), StandardCharsets.UTF_8), '\t', '\0', '\0', "\n");
        for (String col : columns) {
            if(!colNames.add(col)){
                throw new IOException("Unable to add duplicate column name: " + col);
            }
            outputLine[0] = col;
            colNamesWriter.writeNext(outputLine);
        }
        colNamesWriter.close();

        rowNamesWriter = new CSVWriter(new OutputStreamWriter(new GZIPOutputStream(new DigestOutputStream(new FileOutputStream(rowFile), rowFileDigest)), StandardCharsets.UTF_8), '\t', '\0', '\0', "\n");

        blockIndices = new TLongArrayList(1000);


        matrixFileWriter = new CountingOutputStream(new DigestOutputStream(new BufferedOutputStream(new FileOutputStream(matrixFile), 262144), matrixFileDigest));


    }

    /**
     * Write double array to a row. Can be made faster with dedicated
     * implementation instead of converting to DoubleMatrix1D
     *
     * @param rowName The name of the row
     * @param rowData The row data in double[]
     */
    public final void addRow(final String rowName, final double[] rowData) throws IOException {
        if(!rowNames.add(rowName)){
            throw new IOException("Unable to add duplicate row name: " + rowName);
        }

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
            long v = Double.doubleToLongBits(rowData[c]);
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

    public synchronized final void addRow(final String rowName, final DoubleMatrix1D rowData) throws IOException {

        if(!rowNames.add(rowName)){
            throw new IOException("Unable to add duplicate row name: " + rowName);
        }

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
        //System.out.println("rows: " + rowName + " buffer size" + rowBuffer.length + " bytes per row " + bytesPerRow);
        blockCompressionWriter.write(rowBuffer, 0, bytesPerRow);

        if (++rowsInCurrentBlock >= rowsPerBlock) {
            closeBlock();
            //don't start new block because unknown if more data will come
        }

    }

    private void startBlock() throws IOException {
        //System.out.println("Block start: " + matrixFileWriter.getByteCount());
        blockIndices.add(matrixFileWriter.getByteCount());//Set the start byte for this block
        blockCompressionWriter = new FramedLZ4CompressorOutputStream(matrixFileWriter);
        rowsInCurrentBlock = 0;
    }

    private void closeBlock() throws IOException {
        blockCompressionWriter.finish();
        rowsInCurrentBlock = -1;
    }

    public synchronized final void close() throws IOException {
        //here write row indices to end of file
        //Then number of row and columns
        //Finally start of row indices end block

        if (rowsInCurrentBlock > 0) {//>0 is correct this is count not index
            closeBlock();
        }

        final long startOfIndexBlock = matrixFileWriter.getByteCount();

        //System.out.println("Index block start: " + matrixFileWriter.getByteCount());


        final int numberOfBlocks = blockIndices.size();

        //System.out.println("Number of block: " + numberOfBlocks);

        final FramedLZ4CompressorOutputStream rowIndicesCompression = new FramedLZ4CompressorOutputStream(matrixFileWriter);
        final DataOutputStream rowIndicesWriter = new DataOutputStream(rowIndicesCompression);


        for (int r = 0; r < numberOfBlocks; ++r) {
            rowIndicesWriter.writeLong(blockIndices.get(r));
        }
        rowIndicesCompression.finish();

        final DataOutputStream metaDataBlockWriter = new DataOutputStream(matrixFileWriter);

        final long startOfMetaDataBlock = matrixFileWriter.getByteCount();
        //System.out.println("Meta block start: " + matrixFileWriter.getByteCount());

        metaDataBlockWriter.writeUTF(datasetName);
        metaDataBlockWriter.writeUTF(dataOnRows);
        metaDataBlockWriter.writeUTF(dataOnCols);

        metaDataBlockWriter.writeLong(Instant.now().getEpochSecond());
        metaDataBlockWriter.writeInt(numberOfRows);
        metaDataBlockWriter.writeInt(numberOfColumns);
        metaDataBlockWriter.writeLong(startOfIndexBlock);
        metaDataBlockWriter.writeLong(startOfMetaDataBlock);
        //System.out.println("Meta block start: " + rowsPerBlock);

        metaDataBlockWriter.writeInt(rowsPerBlock);
        metaDataBlockWriter.write(0);//reserve
        metaDataBlockWriter.write(0);//reserve
        metaDataBlockWriter.write(0);//reserve
        metaDataBlockWriter.write(0b00000001);//set last bit to 1 to indicate usage lz4 frames for compression
        metaDataBlockWriter.write(MAGIC_BYTES);

        //Also closes matrixFileWriter
        metaDataBlockWriter.close();

        rowNamesWriter.close();

        final BufferedWriter checksumFileWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(checksumFile), StandardCharsets.UTF_8));

        addChecksumToFile(matrixFileDigest.digest(), matrixFile.getName(), checksumFileWriter);
        addChecksumToFile(rowFileDigest.digest(), rowFile.getName(), checksumFileWriter);
        addChecksumToFile(colFileDigest.digest(), colFile.getName(), checksumFileWriter);

        checksumFileWriter.close();

    }

    private void addChecksumToFile(final byte[] hash, final String filename, BufferedWriter checksumFileWriter) throws IOException {

        for (byte b : hash) {
            checksumFileWriter.write(String.format("%02x", b));
        }
        checksumFileWriter.append(' ');
        checksumFileWriter.append(' ');
        checksumFileWriter.write(filename);
        checksumFileWriter.append('\n');

    }

    public static <Comparable> void saveDataset(final String path, final DoubleMatrixDataset<? extends Comparable, ? extends Comparable> dataset) throws FileNotFoundException, IOException {
        saveDataset(path, dataset, "", "", "");
    }

    public static <Comparable> void saveDataset(final String path, final DoubleMatrixDataset<? extends Comparable, ? extends Comparable> dataset, final String datasetName, final String dataOnRows, final String dataOnCols) throws FileNotFoundException, IOException {

        final ArrayList<String> colNames = new ArrayList<>(dataset.columns());

        for (Object c : dataset.getColObjects()) {
            colNames.add(c.toString());
        }

        final DoubleMatrixDatasetRowCompressedWriter writer = new DoubleMatrixDatasetRowCompressedWriter(path, colNames, datasetName, dataOnRows, dataOnCols);

        final ArrayList<? extends Comparable> rowNames = dataset.getRowObjects();

        final int totalRows = dataset.rows();
        for (int r = 0; r < totalRows; ++r) {
            writer.addRow(rowNames.get(r).toString(), dataset.getRow(r));
        }

        writer.close();

    }

}
