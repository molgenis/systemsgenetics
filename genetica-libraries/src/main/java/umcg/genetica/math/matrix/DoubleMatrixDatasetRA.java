/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.matrix;

import java.io.*;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author lude, juha
 */
@Deprecated
public class DoubleMatrixDatasetRA<T, U> extends DoubleMatrixDatasetAC {

    private static final Logger LOGGER = Logger.getLogger(DoubleMatrixDatasetRA.class.getName());
    public int nrColsTotal = 0;
    public Set<T> rowsToInclude = null;
    public Set<U> colsToInclude = null;
    public String fileName = null;
    private RandomAccessFile raf;
    private int[] rowIndexToRawRowIndex;
    private int[] colIndexToRawColIndex;

    public DoubleMatrixDatasetRA(String fileName) throws IOException {
        this(fileName, LoadLabels.LOAD_BOTH);
    }

    public DoubleMatrixDatasetRA(String fileName, LoadLabels ll) throws IOException {
        LOGGER.log(Level.INFO, "Loading dataset: {0}", fileName);
        if (fileName.endsWith(".binary")) {
            try {
                loadExpressionDataInBinaryFormat(fileName, ll);
            } catch (ClassNotFoundException ex) {
                throw new IOException(ex); // for backward compatibility we don't throw ClassNotFoundExceptions :I
            }
        } else {
            throw new IOException("Only .binary format supported for random access! Given file: " + fileName);
        }
    }

    public DoubleMatrixDatasetRA(String fileName, Set<T> rowsToInclude) throws IOException {
        this(fileName, rowsToInclude, LoadLabels.LOAD_BOTH);
    }

    public DoubleMatrixDatasetRA(String fileName, Set<T> rowsToInclude, LoadLabels ll) throws IOException {
        this.rowsToInclude = rowsToInclude;
        LOGGER.log(Level.INFO, "Loading dataset: {0}", fileName);
        if (fileName.endsWith(".binary")) {
            try {
                loadExpressionDataInBinaryFormat(fileName, ll);
            } catch (ClassNotFoundException ex) {
                throw new IOException(ex); // for backward compatibility we don't throw ClassNotFoundExceptions :I
            }
        } else {
            throw new IOException("Only .binary format supported for random access! Given file: " + fileName);
        }
    }

    public DoubleMatrixDatasetRA(String fileName, Set<T> rowsToInclude, Set<U> colsToInclude) throws IOException {
        this(fileName, rowsToInclude, colsToInclude, LoadLabels.LOAD_BOTH);
    }

    public DoubleMatrixDatasetRA(String fileName, Set<T> rowsToInclude, Set<U> colsToInclude, LoadLabels ll) throws IOException {
        this.rowsToInclude = rowsToInclude;
        this.colsToInclude = colsToInclude;
        LOGGER.log(Level.INFO, "Loading dataset: {0}", fileName);
        if (fileName.endsWith(".binary")) {
            try {
                loadExpressionDataInBinaryFormat(fileName, ll);
            } catch (ClassNotFoundException ex) {
                throw new IOException(ex);
            }
        } else {
            throw new IOException("Only .binary format supported for random access! Given file: " + fileName);
        }
    }

    private void loadExpressionDataInBinaryFormat(String fileName, LoadLabels ll) throws IOException, ClassNotFoundException {
        //First load the raw binary data:
        this.fileName = fileName;
        File fileBinary = new File(fileName + ".dat");
        int nrRowsThisBinaryFile = -1;
        int nrColsThisBinaryFile = -1;
        this.raf = new RandomAccessFile(fileBinary, "r");
        raf.seek(0);
        byte[] bytes = new byte[4];
        raf.read(bytes, 0, 4);
        nrRowsThisBinaryFile = byteArrayToInt(bytes);
        raf.read(bytes, 0, 4);
        nrColsThisBinaryFile = byteArrayToInt(bytes);

        if (rowsToInclude == null && colsToInclude == null) {
            //We want to load all the data:
            nrRows = nrRowsThisBinaryFile;
            nrCols = nrColsThisBinaryFile;

            //Now load the row and column identifiers from files
            switch (ll) {
                case LOAD_BOTH:
                    loadRowObjects(fileName, nrRows);
                    loadColumnObjects(fileName, nrCols);
                    break;
                case LOAD_ROWS:
                    loadRowObjects(fileName, nrRows);
                    break;
                case LOAD_COLUMNS:
                    loadColumnObjects(fileName, nrCols);
                    break;
            }

            nrColsTotal = nrCols;

        } else {

            //We want to confine the set of probes and samples to a subset. Deal with this in a different way.
            loadRowObjects(fileName, nrRowsThisBinaryFile);
            loadColumnObjects(fileName, nrColsThisBinaryFile);
            nrColsTotal = nrColsThisBinaryFile;

        }
        recalculateHashMaps();
        LOGGER.log(Level.INFO, "Access to binary file:\t{0}\tok, nrRows:\t{1}\tnrCols:\t{2}", new Object[]{fileName, nrRows, nrCols});
    }

    public static List<Object> getRowObjectsOnly(String fileName) throws FileNotFoundException, IOException, ClassNotFoundException {

        File fileRows = new File(fileName + ".rows.ser");
        List<Object> objects = new ArrayList<Object>();
        BufferedInputStream bis = new BufferedInputStream(new FileInputStream(fileRows));
        ObjectInputStream ois = new ObjectInputStream(bis);

        while (bis.available() > 0) {
            objects.add(ois.readObject());
        }
        bis.close();
        ois.close();

        return objects;
    }

    public static List<Object> getColumnObjectsOnly(String fileName) throws FileNotFoundException, IOException, ClassNotFoundException {

        File fileCols = new File(fileName + ".columns.ser");
        List<Object> objects = new ArrayList<Object>();
        BufferedInputStream bis = new BufferedInputStream(new FileInputStream(fileCols));
        ObjectInputStream ois = new ObjectInputStream(bis);

        while (bis.available() > 0) {
            objects.add(ois.readObject());
        }
        bis.close();
        ois.close();

        return objects;
    }

    private void loadRowObjects(String fileName, int nrRowsThisBinaryFile) throws FileNotFoundException, IOException, ClassNotFoundException {
        int[] rowSubsetIndex = new int[nrRowsThisBinaryFile];
        for (int r = 0; r < rowSubsetIndex.length; r++) {
            rowSubsetIndex[r] = -1;
        }
        File fileRows = new File(fileName + ".rows.ser");
        BufferedInputStream bis = new BufferedInputStream(new FileInputStream(fileRows));
        ObjectInputStream ois = new ObjectInputStream(bis);
        int rowIndex = 0;
        rowObjects = new ArrayList<T>();//rowsToInclude.size());
        HashMap<T, Integer> hashRowsPresentAndRequested = new HashMap<T, Integer>();
        while (bis.available() > 0) {
            T rowObject = (T) ois.readObject();
            if (rowsToInclude == null || rowsToInclude.contains(rowObject)) {
                rowSubsetIndex[rowIndex] = hashRowsPresentAndRequested.size();
                hashRowsPresentAndRequested.put(rowObject, rowIndex);
                rowObjects.add(rowObject);
            }
            rowIndex++;
        }
        bis.close();
        ois.close();

        nrRows = hashRowsPresentAndRequested.size();
        rowIndexToRawRowIndex = new int[nrRows];

        bis = new BufferedInputStream(new FileInputStream(fileRows));
        ois = new ObjectInputStream(bis);
        int rowCounter = 0;
        rowIndex = 0;
        while (bis.available() > 0) {
            T rowObject = (T) ois.readObject();
            if (rowsToInclude == null || rowsToInclude.contains(rowObject)) {
                rowIndexToRawRowIndex[rowCounter] = rowIndex;
                rowCounter++;
            }
            rowIndex++;
        }
        bis.close();
        ois.close();

    }

    private void loadColumnObjects(String fileName, int nrColsThisBinaryFile) throws FileNotFoundException, IOException, ClassNotFoundException {

        int[] colSubsetIndex = new int[nrColsThisBinaryFile];
        for (int c = 0; c < colSubsetIndex.length; c++) {
            colSubsetIndex[c] = -1;
        }
        File fileCols = new File(fileName + ".columns.ser");
        BufferedInputStream bis = new BufferedInputStream(new FileInputStream(fileCols));
        ObjectInputStream ois = new ObjectInputStream(bis);
        int colIndex = 0;
        colObjects = new ArrayList<U>();//colsToInclude.size());
        HashMap<U, Integer> hashColsPresentAndRequested = new HashMap<U, Integer>();
        while (bis.available() > 0) {
            U colObject = (U) ois.readObject();
            if (colsToInclude == null || colsToInclude.contains(colObject)) {
                colSubsetIndex[colIndex] = hashColsPresentAndRequested.size();
                hashColsPresentAndRequested.put(colObject, colIndex);
                colObjects.add(colObject);
            }
            colIndex++;
        }
        bis.close();
        ois.close();

        nrCols = hashColsPresentAndRequested.size();
        colIndexToRawColIndex = new int[nrCols];

        bis = new BufferedInputStream(new FileInputStream(fileCols));
        ois = new ObjectInputStream(bis);
        int colCounter = 0;
        colIndex = 0;
        while (bis.available() > 0) {
            U colObject = (U) ois.readObject();
            if (colsToInclude == null || colsToInclude.contains(colObject)) {
                colIndexToRawColIndex[colCounter] = colIndex;
                colCounter++;
            }
            colIndex++;
        }
        bis.close();
        ois.close();

    }

    public double[] getNextRow() {
        if (nrCols != nrColsTotal) {
            throw new IllegalStateException("Not applicable to datasets with not all columns included. (Columns in file: " + nrColsTotal + ", columns included: " + nrCols + ")");
        }
        double[] values = new double[nrCols];
        int col = 0;
        try {
            for (col = 0; col < nrColsTotal; col++) {
                values[col] = raf.readDouble();
            }
        } catch (IOException e) {
            System.err.println("Can't get element at column " + col + ": " + e.getMessage());
        }
        return values;
    }

    @Override
    public synchronized double[] get(int x) {
        if (nrCols != nrColsTotal) {
            throw new IllegalStateException("Not applicable to datasets with not all columns included. (Columns in file: " + nrColsTotal + ", columns included: " + nrCols + ")");
        }
        int rawRowIndex = x;
        if (rowIndexToRawRowIndex != null) {
            rawRowIndex = rowIndexToRawRowIndex[x];
        }
        double[] values = new double[nrCols];
        int col = 0;
        try {
            raf.seek(8 + (rawRowIndex * (long) nrColsTotal) * 8);
            for (col = 0; col < nrColsTotal; col++) {
                values[col] = raf.readDouble();
            }
        } catch (IOException e) {
            System.err.println("Can't get element (" + rawRowIndex + ", " + col + "): " + e.getMessage());
        }
        return values;
    }

    @Override
    public synchronized double get(int x, int y) {
        int rawRowIndex = x;
        if (rowIndexToRawRowIndex != null) {
            rawRowIndex = rowIndexToRawRowIndex[x];
        }
        int rawColIndex = y;
        if (colIndexToRawColIndex != null) {
            rawColIndex = colIndexToRawColIndex[y];
        }
        double value = Double.NaN;
        try {
            raf.seek(8 + (rawRowIndex * (long) nrColsTotal + rawColIndex) * 8);
            value = raf.readDouble();
        } catch (IOException e) {
            System.err.println("Can't get element (" + rawRowIndex + ", " + rawColIndex + "): " + e.getMessage());
        }
        return value;

    }

    @Override
    public void recalculateHashMaps() {

        if (rowObjects != null) {
            hashRows.clear();
            for (int probeItr = 0; probeItr < nrRows; probeItr++) {
                hashRows.put(rowObjects.get(probeItr), probeItr);
            }
        }
        if (colObjects != null) {
            hashCols.clear();
            for (int sampleItr = 0; sampleItr < nrCols; sampleItr++) {
                hashCols.put(colObjects.get(sampleItr), sampleItr);
            }
        }
    }

    //TODO ByteBuffer
    private byte[] intToByteArray(int value) {
        return new byte[]{(byte) (value >>> 24),
                    (byte) (value >>> 16),
                    (byte) (value >>> 8),
                    (byte) value};
    }

    private int byteArrayToInt(byte[] b) {
        return (b[0] << 24)
                + ((b[1] & 0xff) << 16)
                + ((b[2] & 0xff) << 8)
                + (b[3] & 0xff);
    }
}
