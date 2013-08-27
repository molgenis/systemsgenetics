/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.matrix2;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseLargeDoubleMatrix2D;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author MarcJan
 */
public abstract class DoubleMatrixDataset<R, C> extends DoubleMatrix2D {

    static final Exception doubleMatrixDatasetNonUniqueHeaderException = new Exception("Tried to use a non-unique header set in an identifier HashMap");
    static final Logger LOGGER = Logger.getLogger(DoubleMatrixDataset.class.getName());
    protected LinkedHashMap<R, Integer> hashRows;
    protected LinkedHashMap<C, Integer> hashCols;

    public DoubleMatrixDataset() {
        hashRows = new LinkedHashMap<R, Integer>();
        hashCols = new LinkedHashMap<C, Integer>();
    }

    public DoubleMatrixDataset(LinkedHashMap<R, Integer> hashRows, LinkedHashMap<C, Integer> hashCols) {
        this.hashRows = hashRows;
        this.hashCols = hashCols;
    }

    public abstract DoubleMatrix2D getMatrix();

    public static DoubleMatrixDataset<String, String> loadDoubleData(String fileName) throws IOException, Exception {
        return loadDoubleData(fileName, "\t");
    }

    public static DoubleMatrixDataset<String, String> loadDoubleData(String fileName, String delimiter) throws IOException, Exception {

        Pattern splitPatern = Pattern.compile(delimiter);

        int columnOffset = 1;

        int[] colIndex;
        TextFile in = new TextFile(fileName, TextFile.R);
        String str = in.readLine(); // header
        String[] data = splitPatern.split(str);

        int tmpCols = (data.length - columnOffset);

        LinkedHashMap<String, Integer> colMap = new LinkedHashMap<String, Integer>((int) Math.ceil(tmpCols / 0.75));

        colIndex = new int[tmpCols];
        for (int s = 0; s < tmpCols; s++) {
            String colName = data[s + columnOffset];
            if (!colMap.containsKey(colName)) {
                colMap.put(colName, s);
            } else {
                LOGGER.warning("Duplicated column name!");
                throw (doubleMatrixDatasetNonUniqueHeaderException);
            }
            colIndex[s] = s + columnOffset;
        }

        int tmpRows = 0;

        while (in.readLine() != null) {
            tmpRows++;
        }
        in.close();

        double[][] initialMatrix = new double[tmpRows][tmpCols];

        in.open();
        in.readLine(); // read header
        int row = 0;

        LinkedHashMap<String, Integer> rowMap = new LinkedHashMap<String, Integer>((int) Math.ceil(tmpRows / 0.75));

        boolean correctData = true;
        while ((str = in.readLine()) != null) {
            data = splitPatern.split(str);

            if (!rowMap.containsKey(data[0])) {
                rowMap.put(data[0], row);
            } else {
                LOGGER.warning("Duplicated row name!");
                throw (doubleMatrixDatasetNonUniqueHeaderException);
            }

            for (int s = 0; s < tmpCols; s++) {
                double d;
                try {
                    d = Double.parseDouble(data[s + columnOffset]);
                } catch (NumberFormatException e) {
                    correctData = false;
                    d = Double.NaN;
                }
                initialMatrix[row][s] = d;
            }
            row++;
        }
        if (!correctData) {
            LOGGER.warning("Your data contains NaN/unparseable values!");
        }
        in.close();

        DoubleMatrixDataset<String, String> dataset;

        if ((tmpRows * tmpCols) < (Integer.MAX_VALUE - 2)) {
            dataset = new SmallDoubleMatrixDataset<String, String>(new DenseDoubleMatrix2D(initialMatrix), rowMap, colMap);
        } else {
            DenseLargeDoubleMatrix2D matrix = new DenseLargeDoubleMatrix2D(tmpRows, tmpCols);
            matrix.assign(initialMatrix);
            dataset = new LargeDoubleMatrixDataset<String, String>(matrix, rowMap, colMap);
        }

        LOGGER.log(Level.INFO, "''{0}'' has been loaded, nrRows: {1} nrCols: {2}", new Object[]{fileName, dataset.rows(), dataset.columns()});
        return dataset;
    }

    public static DoubleMatrixDataset<String, String> createMatrixDataset(int rows, int cols) {
        if ((rows * cols) < (Integer.MAX_VALUE - 2)) {
            return new SmallDoubleMatrixDataset<String, String>();
        } else {
            return new LargeDoubleMatrixDataset<String, String>();
        }
    }

    public void saveLowMemory(String fileName) throws IOException {
        TextFile out = new TextFile(fileName, TextFile.W);

        ArrayList<C> colObjects = new ArrayList<C>(hashCols.keySet());
        ArrayList<R> rowObjects = new ArrayList<R>(hashRows.keySet());

        out.append('-');
        for (int s = 0; s < columns; s++) {
            out.append('\t');
            out.append(colObjects.get(s).toString());
        }
        out.append('\n');

        for (int r = 0; r < rows; r++) {
            out.append(rowObjects.get(r).toString());
            DoubleMatrix1D rowInfo = getMatrix().viewRow(r);
            for (int s = 0; s < rowInfo.size(); s++) {
                out.append('\t');
                out.append(String.valueOf(rowInfo.get(s)));
            }
            out.append('\n');
        }
        out.close();
    }

    public void save(String fileName) throws IOException {
        TextFile out = new TextFile(fileName, TextFile.W);

        ArrayList<C> colObjects = new ArrayList<C>(hashCols.keySet());
        ArrayList<R> rowObjects = new ArrayList<R>(hashRows.keySet());

        out.append('-');
        for (int s = 0; s < getMatrix().columns(); s++) {
            out.append('\t');
            out.append(colObjects.get(s).toString());
        }
        out.append('\n');
        double[][] rawData = getMatrix().toArray();

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
        out.close();
    }

    //Getters and setters
    public LinkedHashMap<R, Integer> getHashRows() {
        return hashRows;
    }

    public void setHashRows(LinkedHashMap<R, Integer> hashRows) {
        this.hashRows = hashRows;
    }

    public LinkedHashMap<C, Integer> getHashCols() {
        return hashCols;
    }

    public void setHashCols(LinkedHashMap<C, Integer> hashCols) {
        this.hashCols = hashCols;
    }

    public ArrayList<R> getRowObjects() {
        return new ArrayList<R>(hashRows.keySet());
    }

    public void setRowObjects(ArrayList<R> arrayList) throws Exception {
        LinkedHashMap<R, Integer> newHashRows = new LinkedHashMap<R, Integer>((int) Math.ceil(arrayList.size() / 0.75));
        int i = 0;
        for (R s : arrayList) {
            if (!newHashRows.containsKey(s)) {
                newHashRows.put(s, i);
            } else {
                System.out.println("Error, new row names contains dupilcates.");
                throw (doubleMatrixDatasetNonUniqueHeaderException);
            }
            i++;
        }

        this.hashRows = newHashRows;
    }

    public ArrayList<C> getColObjects() {
        return new ArrayList<C>(hashCols.keySet());
    }

    public void setColObjects(ArrayList<C> arrayList) throws Exception {
        LinkedHashMap<C, Integer> newHashCols = new LinkedHashMap<C, Integer>((int) Math.ceil(arrayList.size() / 0.75));
        int i = 0;
        for (C s : arrayList) {
            if (!newHashCols.containsKey(s)) {
                newHashCols.put(s, i);
            } else {
                System.out.println("Error, new column names contains dupilcates.");
                throw (doubleMatrixDatasetNonUniqueHeaderException);
            }
            i++;
        }
        this.hashCols = newHashCols;
    }

    public abstract void setMatrix(double[][] Matrix);
}
