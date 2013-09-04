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
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.LinkedHashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author MarcJan
 */
public abstract class DoubleMatrixDataset<R extends Comparable, C extends Comparable> extends DoubleMatrix2D {

    static final IOException doubleMatrixDatasetNonUniqueHeaderException = new IOException("Tried to use a non-unique header set in an identifier HashMap");
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

        TextFile in = new TextFile(fileName, TextFile.R);
        String str = in.readLine(); // header
        String[] data = splitPatern.split(str);

        int tmpCols = (data.length - columnOffset);

        LinkedHashMap<String, Integer> colMap = new LinkedHashMap<String, Integer>((int) Math.ceil(tmpCols / 0.75));

        for (int s = 0; s < tmpCols; s++) {
            String colName = data[s + columnOffset];
            if (!colMap.containsKey(colName)) {
                colMap.put(colName, s);
            } else {
                LOGGER.warning("Duplicated column name!");
                throw (doubleMatrixDatasetNonUniqueHeaderException);
            }
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
            } else {
                LOGGER.warning("Duplicated row name!");
                throw (doubleMatrixDatasetNonUniqueHeaderException);
            }           
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
    
    public static DoubleMatrixDataset<String, String> loadSubsetOfDoubleData(String fileName, String delimiter, HashSet<String> desiredRows, HashSet<String> desiredCols) throws IOException, Exception {
        
        LinkedHashSet<Integer> desiredColPos = new LinkedHashSet<Integer>();
        
        Pattern splitPatern = Pattern.compile(delimiter);

        int columnOffset = 1;

        TextFile in = new TextFile(fileName, TextFile.R);
        String str = in.readLine(); // header
        String[] data = splitPatern.split(str);

        int tmpCols = (data.length - columnOffset);

        LinkedHashMap<String, Integer> colMap = new LinkedHashMap<String, Integer>((int) Math.ceil(tmpCols / 0.75));

        int storedCols = 0;
        for (int s = 0; s < tmpCols; s++) {
            String colName = data[s + columnOffset];
            if (!colMap.containsKey(colName) && (desiredCols == null || desiredCols.contains(colName) || desiredCols.isEmpty())) {
                colMap.put(colName, storedCols);
                desiredColPos.add(storedCols);
                storedCols++;
            } else if(colMap.containsKey(colName)){
                LOGGER.warning("Duplicated column name!");
                System.out.println("Tried to add: "+colName);
                throw (doubleMatrixDatasetNonUniqueHeaderException);
            }
        }

        LinkedHashSet<Integer> desiredRowPos = new LinkedHashSet<Integer>();
        int rowsToStore = 0;
        int totalRows = 0;

        while (in.readLine() != null) {
            String[] info = splitPatern.split(str);
            if(desiredRows == null || desiredRows.contains(info[0]) || desiredRows.isEmpty()){
                rowsToStore++;
                desiredRowPos.add(totalRows);
            }
            totalRows++;
        }
        in.close();

        double[][] initialMatrix = new double[rowsToStore][storedCols];

        in.open();
        in.readLine(); // read header
        int storingRow = 0;
        totalRows = 0;
        LinkedHashMap<String, Integer> rowMap = new LinkedHashMap<String, Integer>((int) Math.ceil(rowsToStore / 0.75));

        boolean correctData = true;
        while ((str = in.readLine()) != null) {
            if(desiredRowPos.contains(totalRows)){
                data = splitPatern.split(str);
                if (!rowMap.containsKey(data[0])) {
                    rowMap.put(data[0], storingRow);
                    for (int s : desiredColPos) {
                        double d;
                        try {
                            d = Double.parseDouble(data[s + columnOffset]);
                        } catch (NumberFormatException e) {
                            correctData = false;
                            d = Double.NaN;
                        }
                        initialMatrix[storingRow][s] = d;
                    }
                    storingRow++;
                } else if(rowMap.containsKey(data[0])){
                    LOGGER.warning("Duplicated row name!");
                    System.out.println("Tried to add: "+data[0]);
                    throw (doubleMatrixDatasetNonUniqueHeaderException);
                }
            }
            totalRows++;
        }
        if (!correctData) {
            LOGGER.warning("Your data contains NaN/unparseable values!");
        }
        in.close();

        DoubleMatrixDataset<String, String> dataset;

        if ((rowsToStore * tmpCols) < (Integer.MAX_VALUE - 2)) {
            dataset = new SmallDoubleMatrixDataset<String, String>(new DenseDoubleMatrix2D(initialMatrix), rowMap, colMap);
        } else {
            DenseLargeDoubleMatrix2D matrix = new DenseLargeDoubleMatrix2D(rowsToStore, tmpCols);
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
            out.append(rowObjects.get(p).toString());
            
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
    
    /**
     * Order columns
     *
     * @param dataset DoubleMatrixDataset Expression matrix
     */
    public void OrderOnColumns() {
        LinkedHashMap<C, Integer> newColHash = new LinkedHashMap<C, Integer>((int) Math.ceil(this.columns() / 0.75));
        ArrayList<C> names = this.getColObjects();
        Collections.sort(names);
        
        int pos = 0;
        for(C name : names){
            newColHash.put(name, pos);
            pos++;
        }
        reorderCols(newColHash);
    }
    
    /**
     * Order rows
     *
     * @param dataset DoubleMatrixDataset Expression matrix
     */
    public void OrderOnRows() {
        LinkedHashMap<R, Integer> newRowHash = new LinkedHashMap<R, Integer>((int) Math.ceil(this.rows() / 0.75));
        ArrayList<R> names = this.getRowObjects();
        Collections.sort(names);
        
        int pos = 0;
        for(R name : names){
            newRowHash.put(name, pos);
            pos++;
        }
        reorderRows(newRowHash);

    }

    public abstract void setMatrix(double[][] Matrix);
    
    public abstract void reorderRows(LinkedHashMap<R, Integer> mappingIndex);
    
    public abstract void reorderCols(LinkedHashMap<C, Integer> mappingIndex);
    
    public abstract DoubleMatrixDataset<C,R> viewDice();
    
    //Fixed like in parallel colt.
    @Override
    protected DoubleMatrix1D like1D(int size, int offset, int stride) {
        throw new InternalError(); // should never get called
    }
    
    @Override
    public int hashCode() {
        return 1;
    }
}
