/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.matrix;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map.Entry;
import java.util.StringTokenizer;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author MarcJan
 */
public class DoubleMatrixDataset2<T, U> {

    private static final Logger LOGGER = Logger.getLogger(DoubleMatrixDataset2.class.getName());
    private LinkedHashMap<T, Integer> hashRows = null;
    private LinkedHashMap<U, Integer> hashCols = null;
    private int nrRows;
    private int nrCols;
    private DenseDoubleMatrix2D Matrix = null;
    private Pattern splitPatern;

    // we need this constructor to be able to use the handy save function of this class :D
    public DoubleMatrixDataset2() {
    }

    public DoubleMatrixDataset2(int nrRows, int nrCols) {
        this(nrRows, nrCols, null);
    }

    public DoubleMatrixDataset2(int nrRows, int nrCols, Double initialValue) {

        this.nrRows = nrRows;
        this.nrCols = nrCols;
        // runtime type of the arrays will be Object[] but they can only contain T and U elements
        this.hashRows = new LinkedHashMap<T, Integer>((int) Math.ceil(nrRows / 0.75));
        this.hashCols = new LinkedHashMap<U, Integer>((int) Math.ceil(nrCols / 0.75));
        
            this.Matrix = new DenseDoubleMatrix2D(nrRows, nrCols);
        
        if(initialValue!=null){
            Matrix.assign(initialValue);
        }
    }

    public DoubleMatrixDataset2(String fileName) throws IOException {
        this(fileName, "\t");
    }

    public DoubleMatrixDataset2(String fileName, String ll) throws IOException {
        loadExpressionData(fileName, ll);
    }

    private void loadExpressionData(String fileName, String delimiter) throws IOException {
        splitPatern = Pattern.compile(delimiter);

        int columnOffset = 1;

        int[] colIndex;
        TextFile in = new TextFile(fileName, TextFile.R);
        String str = in.readLine(); // header
        String[] data = splitPatern.split(str);

        nrCols = data.length - columnOffset;

        hashCols = new LinkedHashMap<U, Integer>((int) Math.ceil(nrCols / 0.75));

        colIndex = new int[nrCols];
        for (int s = 0; s < nrCols; s++) {
            String colName = data[s + columnOffset];
            if(!hashCols.containsKey((U)colName)){
                hashCols.put((U)colName, s);
            } else {
                LOGGER.warning("Duplicated column name!");
                System.exit(0);
            }
            colIndex[s] = s + columnOffset;
        }

        nrRows = 0;

        while (in.readLine() != null) {
            nrRows++;
        }
        in.close();

        double[][] initialMatrix = new double[nrRows][nrCols];
        in.open();
        in.readLine(); // read header
        int row = 0;

        hashRows = new LinkedHashMap<T, Integer>((int) Math.ceil(nrRows / 0.75));

        boolean correctData = true;
        while ((str = in.readLine()) != null) {
            data = splitPatern.split(str);
            
            if(!hashRows.containsKey((T)data[0])){
                hashRows.put((T) data[0], row);
            } else {
                LOGGER.warning("Duplicated row name!");
                System.exit(0);
            }

            for (int s = 0; s < nrCols; s++) {
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

        Matrix = new DenseDoubleMatrix2D(initialMatrix);

        LOGGER.log(Level.INFO, "''{0}'' has been loaded, nrRows: {1} nrCols: {2}", new Object[]{fileName, nrRows, nrCols});
    }

    private void loadExpressionDataTokenizer(String fileName, String delimiter) throws IOException {
        splitPatern = Pattern.compile(delimiter);

        int columnOffset = 1;

        int[] colIndex;
        TextFile in = new TextFile(fileName, TextFile.R);
        String str = in.readLine(); // header
        String[] data = splitPatern.split(str);

        nrCols = data.length - columnOffset;

        hashCols = new LinkedHashMap<U, Integer>((int) Math.ceil(nrCols / 0.75));

        colIndex = new int[nrCols];
        for (int s = 0; s < nrCols; s++) {
            String colName = data[s + columnOffset];
            if(!hashCols.containsKey((U) colName)){
                hashCols.put((U) colName, s);
            } else {
                LOGGER.warning("Duplicated column name!");
                System.exit(0);
            }
            colIndex[s] = s + columnOffset;
        }

        nrRows = 0;

        while (in.readLine() != null) {
            nrRows++;
        }
        in.close();

        Matrix = new DenseDoubleMatrix2D(nrRows, nrCols);

        in.open();
        in.readLine(); // read header
        int row = 0;

        hashRows = new LinkedHashMap<T, Integer>((int) Math.ceil(nrRows / 0.75));

        boolean correctData = true;

        while ((str = in.readLine()) != null) {
            StringTokenizer st = new StringTokenizer(str, delimiter);
            int col = 0;
            while (st.hasMoreTokens()) {
                if (col != 0) {
                    double d;
                    try {
                        d = Double.parseDouble(st.nextToken());
                    } catch (NumberFormatException e) {
                        correctData = false;
                        d = Double.NaN;
                    }
                    Matrix.setQuick(row, col, d);
                } else {
                    String key = st.nextToken();
                    if(!hashRows.containsKey((T) key)){
                        hashRows.put((T)key, row);
                    } else {
                        LOGGER.warning("Duplicated row name!");
                        System.exit(0);
                    }
                }
            }
            row++;
        }
        if (!correctData) {
            LOGGER.warning("Your data contains NaN/unparseable values!");
        }
        in.close();
        LOGGER.log(Level.INFO, "''{0}'' has been loaded, nrRows: {1} nrCols: {2}", new Object[]{fileName, nrRows, nrCols});
    }

    //Dit is allemaal al geimporteerd.
    
    public void transposeDataset() {
        this.Matrix = (DenseDoubleMatrix2D) Matrix.viewDice();
        
        nrRows = Matrix.rows();
        nrCols = Matrix.columns();
        
        LinkedHashMap<U, Integer> tmp  = new LinkedHashMap<U, Integer>((int) Math.ceil(nrRows / 0.75));
        for(Entry<U, Integer> entry : hashCols.entrySet()){
            tmp.put((U) entry.getKey(), entry.getValue());
        }
        
        LinkedHashMap<T, Integer> tmp2  = new LinkedHashMap<T, Integer>((int) Math.ceil(nrCols / 0.75));
        for(Entry<T, Integer> entry : hashRows.entrySet()){
            tmp2.put((T) entry.getKey(), entry.getValue());
        }
        
        hashCols = (LinkedHashMap<U, Integer>) tmp2;
        hashRows = (LinkedHashMap<T, Integer>) tmp;
        
    }

    public void saveLowMemory(String fileName) throws IOException {
        TextFile out = new TextFile(fileName, TextFile.W);
        
        ArrayList<String> colObjects = (ArrayList<String>) hashCols.keySet();
        ArrayList<String> rowObjects = (ArrayList<String>) hashRows.keySet();
        
        out.append('-');
        for (int s = 0; s < Matrix.columns(); s++) {
            out.append('\t');
            out.append(colObjects.get(s));
        }
        out.append('\n');

        for (int r = 0; r < Matrix.rows(); r++) {
            out.append(rowObjects.get(r));
            DoubleMatrix1D rowInfo = Matrix.viewRow(r);
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
        
        ArrayList<String> colObjects = new ArrayList(hashCols.keySet());
        ArrayList<String> rowObjects = new ArrayList(hashRows.keySet());
        
        out.append('-');
        for (int s = 0; s < Matrix.columns(); s++) {
            out.append('\t');
            out.append(colObjects.get(s));
        }
        out.append('\n');
        double[][] rawData = Matrix.toArray();

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
    public LinkedHashMap<T, Integer> getHashRows() {
        return hashRows;
    }

    public void setHashRows(LinkedHashMap<T, Integer> hashRows) {
        this.hashRows = hashRows;
    }

    public LinkedHashMap<U, Integer> getHashCols() {
        return hashCols;
    }

    public void setHashCols(LinkedHashMap<U, Integer> hashCols) {
        this.hashCols = hashCols;
    }

    public ArrayList<T> getRowObjects() {
        return new ArrayList<T>(hashRows.keySet());
    }
    
    public void setRowObjects(ArrayList<T> arrayList) {
        LinkedHashMap<T, Integer> newHashRows = new LinkedHashMap<T, Integer>((int) Math.ceil(arrayList.size() / 0.75));
        int i = 0;
        for(T s: arrayList){
            if(!newHashRows.containsKey(s)){
                newHashRows.put(s, i);
            } else {
                System.out.println("Error, new row names contains dupilcates.");
                System.exit(-1);
            }
            i++;
        }
        
        this.hashRows = newHashRows;
    }

    public ArrayList<U> getColObjects() {
        return new ArrayList<U>(hashCols.keySet());
    }
    
    public void setColObjects(ArrayList<U> arrayList) {
        LinkedHashMap<U, Integer> newHashCols = new LinkedHashMap<U, Integer>((int) Math.ceil(arrayList.size() / 0.75));
        int i = 0;
        for(U s: arrayList){
            if(!newHashCols.containsKey(s)){
                newHashCols.put(s, i);
            } else {
                System.out.println("Error, new row names contains dupilcates.");
                System.exit(-1);
            }
            i++;
        }
         this.hashCols = newHashCols;
    }

    public int getNrRows() {
        return nrRows;
    }

    public void setNrRows(int nrRows) {
        this.nrRows = nrRows;
    }

    public int getNrCols() {
        return nrCols;
    }

    public void setNrCols(int nrCols) {
        this.nrCols = nrCols;
    }

    public DenseDoubleMatrix2D getMatrix() {
        return Matrix;
    }

    public void setMatrix(DenseDoubleMatrix2D Matrix) {
        this.Matrix = Matrix;
    }
    
    public void setMatrix(double[][] Matrix) {
        this.Matrix = new DenseDoubleMatrix2D(Matrix);
    }
}
