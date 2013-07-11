/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.matrix;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.concurrent.DoubleParseTask;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDatasetAC.LoadLabels;

/**
 *
 * @author lude, juha
 */
public class DoubleMatrixDataset<T, U> extends DoubleMatrixDatasetAC<T, U> {
//    public enum LoadLabels {
//
//        LOAD_BOTH, LOAD_ROWS, LOAD_COLUMNS, DONT_LOAD
//    };

    private static final Logger LOGGER = Logger.getLogger(DoubleMatrixDataset.class.getName());
    public double[][] rawData = null;
    //public int nrRows = 0;
    //public int nrCols = 0;
    //public List<T> rowObjects = null;
    //public List<U> colObjects = null;
    public Set<T> rowsToInclude = null;
    public Set<U> colsToInclude = null;
    //public Map<T, Integer> hashRows = new HashMap<T, Integer>();
    //public Map<U, Integer> hashCols = new HashMap<U, Integer>();
    public String fileName = null;
    private int columnOffset = -1;

    // we need this constructor to be able to use the handy save function of this class :D
    public DoubleMatrixDataset() {
    }

    public DoubleMatrixDataset(String fileName) throws IOException {
        this(fileName, LoadLabels.LOAD_BOTH);
    }

    public DoubleMatrixDataset(String fileName, LoadLabels ll) throws IOException {
        LOGGER.log(Level.INFO, "Loading dataset: {0}", fileName);
        if (fileName.endsWith(".binary")) {
            try {
                loadExpressionDataInBinaryFormat(fileName, ll);
            } catch (ClassNotFoundException ex) {
                throw new IOException(ex); // for backward compatibility we don't throw ClassNotFoundExceptions :I
            }
        } else {
            loadExpressionData(fileName, "\t");
        }
    }

    public DoubleMatrixDataset(String fileName, Set<T> rowsToInclude) throws IOException {
        this(fileName, rowsToInclude, LoadLabels.LOAD_BOTH);
    }

    public DoubleMatrixDataset(String fileName, Set<T> rowsToInclude, LoadLabels ll) throws IOException {
        this.rowsToInclude = rowsToInclude;
        LOGGER.log(Level.INFO, "Loading dataset: {0}", fileName);
        if (fileName.endsWith(".binary")) {
            try {
                loadExpressionDataInBinaryFormat(fileName, ll);
            } catch (ClassNotFoundException ex) {
                throw new IOException(ex); // for backward compatibility we don't throw ClassNotFoundExceptions :I
            }
        } else {
            loadExpressionData(fileName, "\t");
        }
    }

    public DoubleMatrixDataset(String fileName, Set<T> rowsToInclude, Set<U> colsToInclude) throws IOException {
        this(fileName, rowsToInclude, colsToInclude, LoadLabels.LOAD_BOTH);
    }

    public DoubleMatrixDataset(String fileName, Set<T> rowsToInclude, Set<U> colsToInclude, int columnOffset) throws IOException {
        this(fileName, rowsToInclude, colsToInclude, LoadLabels.LOAD_BOTH, columnOffset);
    }

    public DoubleMatrixDataset(String fileName, Set<T> rowsToInclude, Set<U> colsToInclude, LoadLabels ll) throws IOException {
        this(fileName, rowsToInclude, colsToInclude, LoadLabels.LOAD_BOTH, 1);
    }

    public DoubleMatrixDataset(String fileName, Set<T> rowsToInclude, Set<U> colsToInclude, LoadLabels ll, int columnOffset) throws IOException {
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
            loadExpressionData(fileName, "\t", columnOffset);
        }
    }

    // TODO txt file gson
    public DoubleMatrixDataset(String fileName, String delimiter) throws IOException {
        if (!fileName.endsWith(".txt")) {
            throw new IllegalArgumentException("File type must be .txt when delimiter is given (given filename: " + fileName + ")");
        }
        LOGGER.log(Level.INFO, "Loading dataset: {0}", fileName);
        loadExpressionData(fileName, delimiter);
    }

    public DoubleMatrixDataset(double[][] data) {
        this(data.length, data[0].length);
        this.rawData = data;
    }

    public DoubleMatrixDataset(double[][] data, List<U> rowNames, List<T> columnNames) {
        if (data == null || rowNames == null || columnNames == null) {
            throw new IllegalArgumentException("Invalid constructor invocation: supply raw data, column and row names");
        }
        if (data.length != rowNames.size()) {
            throw new IllegalArgumentException("Length of matrix does not correspond to number of row names: " + data.length + " vs " + rowNames.size());
        }
        if (data[0].length != columnNames.size()) {
            throw new IllegalArgumentException("Length of matrix does not correspond to number of col names: " + data[0].length + " vs " + columnNames.size());
        }
        rawData = data;
        rowObjects = (List<T>) rowNames;
        colObjects = (List<U>) columnNames;
        nrCols = colObjects.size();
        nrRows = rowObjects.size();
        recalculateHashMaps();
    }

    public DoubleMatrixDataset(int nrRows, int nrCols) {
        this.nrRows = nrRows;
        this.nrCols = nrCols;
        // runtime type of the arrays will be Object[] but they can only contain T and U elements
        this.rowObjects = new ArrayList<T>(nrRows);
        for (int i = 0; i < nrRows; i++) {
            this.rowObjects.add(null);
        }
        this.colObjects = new ArrayList<U>(nrCols);
        for (int i = 0; i < nrCols; i++) {
            this.colObjects.add(null);
        }
        this.rawData = new double[nrRows][nrCols];
    }

    private void loadExpressionDataInBinaryFormat(String fileName, LoadLabels ll) throws IOException, ClassNotFoundException {
        //First load the raw binary data:
        this.fileName = fileName;
        File fileBinary = new File(fileName + ".dat");
        BufferedInputStream in = null;
        int nrRowsThisBinaryFile = -1;
        int nrColsThisBinaryFile = -1;
        in = new BufferedInputStream(new FileInputStream(fileBinary));
        byte[] bytes = new byte[4];
        in.read(bytes, 0, 4);
        nrRowsThisBinaryFile = byteArrayToInt(bytes);
        in.read(bytes, 0, 4);
        nrColsThisBinaryFile = byteArrayToInt(bytes);

        if (rowsToInclude == null && colsToInclude == null) {
            //We want to load all the data:
            nrRows = nrRowsThisBinaryFile;
            nrCols = nrColsThisBinaryFile;
            rawData = new double[nrRows][nrCols];

            //Now load the row and column identifiers from files
            switch (ll) {
                case LOAD_BOTH:
                    loadRowObjects(fileName);
                    loadColumnObjects(fileName);
                    break;
                case LOAD_ROWS:
                    loadRowObjects(fileName);
                    break;
                case LOAD_COLUMNS:
                    loadColumnObjects(fileName);
                    break;
            }

            byte[] buffer = new byte[nrCols * 8];
            long bits = 0;
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

                    rawData[row][col] = Double.longBitsToDouble(bits);
                    bufferLoc += 8;
                }
            }
            in.close();
        } else {

            //We want to confine the set of probes and samples to a subset. Deal with this in a different way.
            int[] rowSubsetIndex = loadRowObjects(fileName, nrRowsThisBinaryFile);
            int[] colSubsetIndex = loadColumnObjects(fileName, nrColsThisBinaryFile);

            //Now load the binary data:
            rawData = new double[nrRows][nrCols];
            byte[] buffer = new byte[nrColsThisBinaryFile * 8];
            long bits = 0;
            for (int row = 0; row < nrRowsThisBinaryFile; row++) {
                in.read(buffer, 0, nrColsThisBinaryFile * 8);
                int bufferLoc = 0;
                for (int col = 0; col < nrColsThisBinaryFile; col++) {
                    bits = (long) (0xff & buffer[bufferLoc + 7])
                            | (long) (0xff & buffer[bufferLoc + 6]) << 8
                            | (long) (0xff & buffer[bufferLoc + 5]) << 16
                            | (long) (0xff & buffer[bufferLoc + 4]) << 24
                            | (long) (0xff & buffer[bufferLoc + 3]) << 32
                            | (long) (0xff & buffer[bufferLoc + 2]) << 40
                            | (long) (0xff & buffer[bufferLoc + 1]) << 48
                            | (long) (buffer[bufferLoc]) << 56;

                    int rowI = rowSubsetIndex[row];
                    int colI = colSubsetIndex[col];
                    if (rowI != -1 && colI != -1) {
                        rawData[rowI][colI] = Double.longBitsToDouble(bits);
                    }
                    bufferLoc += 8;
                }
            }
            in.close();
        }
        recalculateHashMaps();
        LOGGER.log(Level.INFO, "Binary file ''{0}'' has been loaded, nrRows: {1} nrCols: {2}", new Object[]{fileName, nrRows, nrCols});
    }

    @SuppressWarnings("element-type-mismatch") // generics don't work too well here with text files if the row/col id's are just plain text (use gson or something to represent objects if needed)
    public void loadExpressionData(String fileName, String delimiter, int columnOffset) throws IOException {
        this.columnOffset = columnOffset;
        loadExpressionData(fileName, delimiter);
    }

    @SuppressWarnings("element-type-mismatch") // generics don't work too well here with text files if the row/col id's are just plain text (use gson or something to represent objects if needed)
    private void loadExpressionData(String fileName, String delimiter) throws IOException {
        this.fileName = fileName;

        columnOffset = 1;
        int[] colIndex = null;
        TextFile in = new TextFile(fileName, TextFile.R);
        String str = in.readLine(); // header
        String[] data = str.split(delimiter); // split header


        if (colsToInclude == null) {
            nrCols = data.length - columnOffset;
            colObjects = new ArrayList<U>(nrCols);
            colIndex = new int[nrCols];
            for (int s = 0; s < nrCols; s++) {
                String colName = data[s + columnOffset];
                colObjects.add((U) colName);
                hashCols.put((U) colName, s);
                colIndex[s] = s + columnOffset;
            }
        } else {
            int nrMaxCols = data.length - columnOffset;
            int colCtr = 0;
            for (int s = 0; s < nrMaxCols; s++) {
                String colName = data[s + columnOffset];
                if (colsToInclude.contains(colName)) {
                    colCtr++;
                }
            }

            colObjects = new ArrayList<U>(colCtr);
            colIndex = new int[colCtr]; // this will map from position in memory to position in file...
            nrCols = colCtr;
            colCtr = 0;
            for (int s = 0; s < nrMaxCols; s++) {
                String colName = data[s + columnOffset];
                if (colsToInclude.contains(colName)) {
                    colObjects.add((U) colName);
                    hashCols.put((U) colName, colCtr);  
                    colIndex[colCtr] = s + columnOffset;
                    colCtr++;
                } 
            }

        }



        nrRows = 0;

        ArrayList<String> rowNames = new ArrayList<String>();
        while ((str = in.readLine()) != null) {
            String[] split = str.split(delimiter);
            if (rowsToInclude == null) {
                rowNames.add(split[0]);
                nrRows++;
            } else {
                if (rowsToInclude.contains(split[0])) {
                    rowNames.add(split[0]);
                    nrRows++;
                }
            }
        }
        in.close();

        rowObjects = new ArrayList<T>(nrRows);
        for (int i = 0; i < rowNames.size(); i++) {
            rowObjects.add((T) rowNames.get(i));
        }

        rawData = new double[nrRows][nrCols];

        in.open();
        str = in.readLine(); // read header

        int nrprocs = Runtime.getRuntime().availableProcessors();
//        System.out.println("Using " + nrprocs + " threads for parsing the file");
        ExecutorService threadPool = Executors.newFixedThreadPool(nrprocs);
        CompletionService<Triple<Integer, String, double[]>> pool = new ExecutorCompletionService<Triple<Integer, String, double[]>>(threadPool);

        boolean correctData = true;
        int tasksSubmitted = 0;
        int returnedResults = 0;
        int lnctr = 0;

        Pattern p = Pattern.compile(delimiter);

//        str = in.readLine(); // first header
        while ((str = in.readLine()) != null) {
            DoubleParseTask task = new DoubleParseTask(str, columnOffset, lnctr, colIndex, rowsToInclude, p);
            pool.submit(task);

            tasksSubmitted++;

            if (lnctr % (nrprocs * 2) == 0) {
                while (returnedResults < tasksSubmitted) {
                    try {
                        Triple<Integer, String, double[]> result = pool.take().get();
                        if (result != null) {
                            int rownr = result.getLeft(); //  < 0 when row is not to be included because of hashProbesToInclude.
                            if (rownr >= 0) {
                                double[] doubles = result.getRight();
                                rawData[rownr] = doubles;
                            }
                            result = null;
                            returnedResults++;
                        }
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                }
            }


//            data = str.split(delimiter);
//            if (rowsToInclude == null || rowsToInclude.contains(data[0])) {
//                rowObjects.add((T) new String(data[0].getBytes()));
//                hashRows.put(rowObjects.get(row), row);
//                for (int s = 0; s < nrCols; s++) {
//                    String cell = data[sampleIndex[s] + sampleOffset];
//                    Double d = Double.NaN;
//                    try {
//                        d = Double.parseDouble(cell);
//                    } catch (NumberFormatException e) {
//                        correctData = false;
//                    }
//                    rawData[row][s] = d;
//                }
//                row++;
//            }
//            if (row % 1000 == 0) {
//                System.out.println(row);
//            }
            lnctr++;
        }

        while (returnedResults < tasksSubmitted) {
            try {
                Triple<Integer, String, double[]> result = pool.take().get();
                if (result != null) {
                    int rownr = result.getLeft(); //  < 0 when row is not to be included because of hashProbesToInclude.
                    if (rownr >= 0) {
                        double[] doubles = result.getRight();
                        rawData[rownr] = doubles;
                    }
                    result = null;
                    returnedResults++;
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        threadPool.shutdown();

        in.close();
        recalculateHashMaps();
        LOGGER.log(Level.INFO, "''{0}'' has been loaded, nrRows: {1} nrCols: {2}", new Object[]{fileName, nrRows, nrCols});
    }

    private int[] loadRowObjects(String fileName, int nrRowsThisBinaryFile) throws FileNotFoundException, IOException, ClassNotFoundException {
        int[] rowSubsetIndex = new int[nrRowsThisBinaryFile];
        for (int r = 0; r < rowSubsetIndex.length; r++) {
            rowSubsetIndex[r] = -1;
        }
        File fileRows = new File(fileName + ".rows.ser");
        BufferedInputStream bis = new BufferedInputStream(new FileInputStream(fileRows));
        ObjectInputStream ois = new ObjectInputStream(bis);
        int rowIndex = 0;
        rowObjects = new ArrayList<T>();
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
        return rowSubsetIndex;
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

    private void loadRowObjects(String fileName) throws FileNotFoundException, IOException {

        rowObjects = new ArrayList<T>(nrRows);
        File fileRows = new File(fileName + ".rows.ser");
//            File fileRows = new File(fileName + ".rows.gson");

        BufferedInputStream bis = new BufferedInputStream(new FileInputStream(fileRows));
        ObjectInputStream ois = new ObjectInputStream(bis);

//            TextFile tf = new TextFile(fileName + ".rows.gson", false);
//            String line;
//            Type type = new TypeToken<T>() {
//            }.getType(); // workaround for generics
//            Gson gson = new Gson();
//            while ((line = tf.readLine()) != null) {
//                T fromJson = gson.fromJson(line, type);
//                rowObjects.add(fromJson);
//            }
//            tf.close();

        try {
            while (bis.available() > 0) {
                rowObjects.add((T) ois.readObject());
            }
        } catch (ClassNotFoundException ex) {
            Logger.getLogger(DoubleMatrixDataset.class.getName()).log(Level.SEVERE, "Error with row objects", ex);
            ex.printStackTrace();
        }
        bis.close();
        ois.close();
        if (rowObjects.size() != nrRows) {
            throw new IOException("The number of row objects in " + fileRows.getName() + " doesn't match the number of rows in " + fileName + ".dat");
        }

    }

    private int[] loadColumnObjects(String fileName, int nrColsThisBinaryFile) throws FileNotFoundException, IOException, ClassNotFoundException {

        int[] colSubsetIndex = new int[nrColsThisBinaryFile];
        for (int c = 0; c < colSubsetIndex.length; c++) {
            colSubsetIndex[c] = -1;
        }
        File fileCols = new File(fileName + ".columns.ser");
        BufferedInputStream bis = new BufferedInputStream(new FileInputStream(fileCols));
        ObjectInputStream ois = new ObjectInputStream(bis);
        int colIndex = 0;
        colObjects = new ArrayList<U>();
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
        return colSubsetIndex;
    }

    private void loadColumnObjects(String fileName) throws FileNotFoundException, IOException {

        colObjects = new ArrayList<U>(nrCols);
        File fileCols = new File(fileName + ".columns.ser");
//            File fileCols = new File(fileName + ".columns.gson");

        BufferedInputStream bis = new BufferedInputStream(new FileInputStream(fileCols));
        ObjectInputStream ois = new ObjectInputStream(bis);

//            TextFile tf = new TextFile(fileName + ".rows.gson", false);
//            String line;
//            Type type = new TypeToken<T>() {
//            }.getType(); // workaround for generics
//            Gson gson = new Gson();
//            while ((line = tf.readLine()) != null) {
//                T fromJson = gson.fromJson(line, type);
//                rowObjects.add(fromJson);
//            }
//            tf.close();

        try {
            while (bis.available() > 0) {
                colObjects.add((U) ois.readObject());
            }
        } catch (ClassNotFoundException ex) {
            Logger.getLogger(DoubleMatrixDataset.class.getName()).log(Level.SEVERE, "Error with column objects", ex);
            ex.printStackTrace();
        }
        bis.close();
        ois.close();
        if (colObjects.size() != nrCols) {
            throw new IOException("The number of column objects in " + fileCols.getName() + " doesn't match the number of columns in " + fileName + ".dat");
        }

    }

    public void recalculateHashMaps() {

        if (rowObjects != null) {

            nrRows = rowObjects.size();
            hashRows.clear();
            for (int probeItr = 0; probeItr < nrRows; probeItr++) {
                hashRows.put(rowObjects.get(probeItr), probeItr);
            }
        }
        if (colObjects != null) {
            nrCols = colObjects.size();
            hashCols.clear();
            for (int sampleItr = 0; sampleItr < nrCols; sampleItr++) {
                hashCols.put(colObjects.get(sampleItr), sampleItr);
            }
        }
    }

    public void transposeDataset() {
        rawData = getRawDataTransposed();

        int nrColsTemp = rowObjects.size();
        List<T> colObjsTemp = rowObjects;

        nrRows = colObjects.size();
        rowObjects = new ArrayList<T>();
        for (int i = 0; i < nrRows; i++) {
            rowObjects.add((T) colObjects.get(i));
        }

        nrCols = nrColsTemp;
        colObjects = new ArrayList<U>();

        for (int i = 0; i < nrCols; i++) {
            colObjects.add((U) colObjsTemp.get(i));
        }

        recalculateHashMaps();

    }

    public double[][] getRawDataTransposed() {
        double[][] rawDataTransposed = new double[nrCols][nrRows];
        for (int s = 0; s < nrCols; s++) {
            for (int p = 0; p < nrRows; p++) {
                rawDataTransposed[s][p] = rawData[p][s];
            }
        }
        return rawDataTransposed;
    }

    public void calculateProbeCoexpression(String outputFile) throws IOException {

        standardNormalizeData();

        System.gc();
        System.gc();
        System.gc();
        System.gc();

        SymmetricByteDistanceMatrix matrix = new SymmetricByteDistanceMatrix(nrRows);
        int sampleCountMinusOne = nrCols - 1;
        long[] corrDist = new long[201];
        for (int f = 0; f < nrRows; f++) {
            for (int g = f + 1; g < nrRows; g++) {
                //Calculate correlation:
                double covarianceInterim = 0;
                for (int s = 0; s < nrCols; s++) {
                    covarianceInterim += rawData[f][s] * rawData[g][s];
                }
                double covariance = covarianceInterim / (double) (sampleCountMinusOne);
                double correlation = covariance;
                int corrIndex = (int) Math.round((correlation + 1.0d) * 100.0d);
                corrDist[corrIndex]++;
                matrix.set(f, g, corrIndex);
            }
            if (f % 1000 == 999) {
                LOGGER.log(Level.INFO, "{0} Probes processed", (f + 1));
            }
        }
        matrix.save(new File(outputFile));

        TextFile out = new TextFile(outputFile + "-MatrixProbes.txt", TextFile.W);
        for (int f = 0; f < nrRows; f++) {
            out.write(rowObjects.get(f) + "\t" + rowObjects.get(f) + "\n");
        }
        out.close();
    }

    public void standardNormalizeData() {

        /*
         * System.out.println("\nNormalizing data:"); //Calculate the average
         * expression, when per sample all raw expression levels have been
         * ordered: double[] rankedMean = new double[nrProbes]; for (int
         * probeID=0; probeID<nrProbes; probeID++) { double quantile = ((double)
         * probeID + 1.0d) / ((double) nrProbes + 1d); rankedMean[probeID] =
         * cern.jet.stat.Probability.normalInverse(quantile); }
         *
         * //Iterate through each sample: hgea.math.RankDoubleArray
         * rankDoubleArray = new hgea.math.RankDoubleArray(); for (int s=0;
         * s<nrSamples; s++) { double[] probes = new double[nrProbes]; for (int
         * p=0; p<nrProbes; p++) probes[p]=rawData[p][s]; double[] probesRanked
         * = rankDoubleArray.rank(probes); double[] probesQuantileNormalized =
         * new double[nrProbes]; for (int p=0; p<nrProbes; p++) {
         * probesQuantileNormalized[p] = rankedMean[(int) probesRanked[p]]; }
         * for (int p=0; p<nrProbes; p++) rawData[p][s] = (float)
         * probesQuantileNormalized[p]; }
         */


        LOGGER.info("Setting probe mean to zero and stdev to one for every probe:");
        for (int probeID = 0; probeID < nrRows; probeID++) {
            double vals[] = new double[nrCols];
            System.arraycopy(rawData[probeID], 0, vals, 0, nrCols);
            double mean = JSci.maths.ArrayMath.mean(vals);
            for (int s = 0; s < nrCols; s++) {
                vals[s] -= (double) mean;
            }
            double standardDeviation = JSci.maths.ArrayMath.standardDeviation(vals);
            for (int s = 0; s < nrCols; s++) {
                rawData[probeID][s] = (float) (vals[s] / standardDeviation);
            }
        }

    }

    public void save(String fileName) throws IOException {
        if (fileName.endsWith(".binary") || fileName.endsWith(".dat")) {

            //Create binary file:
            BufferedOutputStream out = null;
            File fileBinary = new File(fileName + ".dat");
            out = new BufferedOutputStream(new FileOutputStream(fileBinary));
            out.write(intToByteArray(nrRows));
            out.write(intToByteArray(nrCols));
            byte[] buffer = new byte[rawData[0].length * 8];
            for (int row = 0; row < rawData.length; row++) { // rows
                int bufferLoc = 0;
                for (int col = 0; col < rawData[0].length; col++) { // columns
                    long bits = Double.doubleToLongBits(rawData[row][col]);
                    buffer[bufferLoc] = (byte) (bits >> 56);
                    buffer[bufferLoc + 1] = (byte) (bits >> 48 & 0xff);
                    buffer[bufferLoc + 2] = (byte) (bits >> 40 & 0xff);
                    buffer[bufferLoc + 3] = (byte) (bits >> 32 & 0xff);
                    buffer[bufferLoc + 4] = (byte) (bits >> 24 & 0xff);
                    buffer[bufferLoc + 5] = (byte) (bits >> 16 & 0xff);
                    buffer[bufferLoc + 6] = (byte) (bits >> 8 & 0xff);
                    buffer[bufferLoc + 7] = (byte) (bits & 0xff);
                    bufferLoc += 8;
                }
                out.write(buffer);
            }
            out.close();
            File fileRows = new File(fileName + ".rows.ser");
            ObjectOutputStream outRows = new ObjectOutputStream(new FileOutputStream(fileRows));
            for (int p = 0; p < rawData.length; p++) {
                outRows.writeObject(rowObjects.get(p));
            }
            outRows.close();
            File fileCols = new File(fileName + ".columns.ser");
            ObjectOutputStream outCols = new ObjectOutputStream(new FileOutputStream(fileCols));
            for (int p = 0; p < rawData[0].length; p++) {
                outCols.writeObject(colObjects.get(p));
            }
            outCols.close();

        } else {
            TextFile out = new TextFile(fileName, TextFile.W);
            out.append('-');
            for (int s = 0; s < rawData[0].length; s++) {
                if (colObjects == null) {
                } else {
                    out.append('\t');
                    out.append(String.valueOf(colObjects.get(s)));
                }

            }
            out.append('\n');

            for (int p = 0; p < rawData.length; p++) {
                if (rowObjects == null) {
                    out.append("");
                } else {
                    out.append(rowObjects.get(p).toString());
                }
                for (int s = 0; s < rawData[p].length; s++) {
                    out.append('\t');
                    out.append(String.valueOf(rawData[p][s]));
                }
                out.append('\n');
            }
            out.close();
        }
    }

    @Override
    public double[] get(int x) {
        return rawData[x];
    }

    @Override
    public double get(int x, int y) {
        return rawData[x][y];
    }

    /**
     *
     * @return transposed dataset with references to the SAME row and column
     * objects and maps, no cloning
     */
    public DoubleMatrixDataset<U, T> getTransposedDataset() {
        DoubleMatrixDataset<U, T> transposed = new DoubleMatrixDataset<U, T>(getRawDataTransposed());
        transposed.rowObjects = colObjects;
        transposed.hashRows = hashCols;
        transposed.colObjects = rowObjects;
        transposed.hashCols = hashRows;
        return transposed;

    }

    public void removeColumnsWithNaNs() {

        Set<Integer> colsWithNaNs = new HashSet<Integer>();
        for (int i = 0; i < nrRows; i++) {
            for (int j = 0; j < nrCols; j++) {
                if (Double.isNaN(rawData[i][j])) {
                    colsWithNaNs.add(j);
                }
            }
        }
        if (!colsWithNaNs.isEmpty()) {
            double[][] newData = new double[nrRows][nrCols - colsWithNaNs.size()];
            for (int i = 0; i < nrRows; i++) {
                int colI = 0;
                for (int j = 0; j < nrCols; j++) {
                    if (!colsWithNaNs.contains(j)) {
                        newData[i][colI] = rawData[i][j];
                        colI++;
                    }
                }
            }
            rawData = newData;
            colObjects.removeAll(colsWithNaNs);
            nrCols -= colsWithNaNs.size();
            recalculateHashMaps();
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

    public double[][] getRawData() {
        return rawData;
    }

    public void setRawData(double[][] rawData) {
        this.rawData = rawData;
    }
}
