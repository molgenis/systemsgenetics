/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.legacy;

//import cern.colt.matrix.tdouble.DoubleMatrix1D;
//import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
//import java.io.IOException;
//import java.util.ArrayList;
//import java.util.LinkedHashMap;
//import java.util.Map.Entry;
//import java.util.StringTokenizer;
//import java.util.logging.Level;
//import java.util.logging.Logger;
//import java.util.regex.Pattern;
//import umcg.genetica.io.text.TextFile;

/**
 *
 * @author MarcJan
 */
public class DoubleMatrixDataset2<T, U> {

//    private static final Logger LOGGER = Logger.getLogger(DoubleMatrixDataset2.class.getName());
//    private LinkedHashMap<T, Integer> hashRows = null;
//    private LinkedHashMap<U, Integer> hashCols = null;
//    private int nrRows;
//    private int nrCols;
//    private DenseDoubleMatrix2D Matrix = null;
//    private Pattern splitPatern;
//
//    // we need this constructor to be able to use the handy save function of this class :D
//    public DoubleMatrixDataset2() {
//    }
//
//    public DoubleMatrixDataset2(int nrRows, int nrCols) {
//        this(nrRows, nrCols, null);
//    }
//
//    public DoubleMatrixDataset2(int nrRows, int nrCols, Double initialValue) {
//
//        this.nrRows = nrRows;
//        this.nrCols = nrCols;
//        // runtime type of the arrays will be Object[] but they can only contain T and U elements
//        this.hashRows = new LinkedHashMap<T, Integer>((int) Math.ceil(nrRows / 0.75));
//        this.hashCols = new LinkedHashMap<U, Integer>((int) Math.ceil(nrCols / 0.75));
//        
//            this.Matrix = new DenseDoubleMatrix2D(nrRows, nrCols);
//        
//        if(initialValue!=null){
//            Matrix.assign(initialValue);
//        }
//    }
//
//    
//    //Dit is allemaal al geimporteerd.
//    
//    public void transposeDataset() {
//        this.Matrix = (DenseDoubleMatrix2D) Matrix.viewDice();
//        
//        nrRows = Matrix.rows();
//        nrCols = Matrix.columns();
//        
//        LinkedHashMap<U, Integer> tmp  = new LinkedHashMap<U, Integer>((int) Math.ceil(nrRows / 0.75));
//        for(Entry<U, Integer> entry : hashCols.entrySet()){
//            tmp.put((U) entry.getKey(), entry.getValue());
//        }
//        
//        LinkedHashMap<T, Integer> tmp2  = new LinkedHashMap<T, Integer>((int) Math.ceil(nrCols / 0.75));
//        for(Entry<T, Integer> entry : hashRows.entrySet()){
//            tmp2.put((T) entry.getKey(), entry.getValue());
//        }
//        
//        hashCols = (LinkedHashMap<U, Integer>) tmp2;
//        hashRows = (LinkedHashMap<T, Integer>) tmp;
//        
//    }
//
//    public void saveLowMemory(String fileName) throws IOException {
//        TextFile out = new TextFile(fileName, TextFile.W);
//        
//        ArrayList<String> colObjects = (ArrayList<String>) hashCols.keySet();
//        ArrayList<String> rowObjects = (ArrayList<String>) hashRows.keySet();
//        
//        out.append('-');
//        for (int s = 0; s < Matrix.columns(); s++) {
//            out.append('\t');
//            out.append(colObjects.get(s));
//        }
//        out.append('\n');
//
//        for (int r = 0; r < Matrix.rows(); r++) {
//            out.append(rowObjects.get(r));
//            DoubleMatrix1D rowInfo = Matrix.viewRow(r);
//            for (int s = 0; s < rowInfo.size(); s++) {
//                out.append('\t');
//                out.append(String.valueOf(rowInfo.get(s)));
//            }
//            out.append('\n');
//        }
//        out.close();
//    }
//
//    public void save(String fileName) throws IOException {
//        TextFile out = new TextFile(fileName, TextFile.W);
//        
//        ArrayList<String> colObjects = new ArrayList(hashCols.keySet());
//        ArrayList<String> rowObjects = new ArrayList(hashRows.keySet());
//        
//        out.append('-');
//        for (int s = 0; s < Matrix.columns(); s++) {
//            out.append('\t');
//            out.append(colObjects.get(s));
//        }
//        out.append('\n');
//        double[][] rawData = Matrix.toArray();
//
//        for (int p = 0; p < rawData.length; p++) {
//            if (rowObjects == null) {
//                out.append("");
//            } else {
//                out.append(rowObjects.get(p).toString());
//            }
//            for (int s = 0; s < rawData[p].length; s++) {
//                out.append('\t');
//                out.append(String.valueOf(rawData[p][s]));
//            }
//            out.append('\n');
//        }
//        out.close();
//        out.close();
//    }
//    
//    //Getters and setters
//    public LinkedHashMap<T, Integer> getHashRows() {
//        return hashRows;
//    }
//
//    public void setHashRows(LinkedHashMap<T, Integer> hashRows) {
//        this.hashRows = hashRows;
//    }
//
//    public LinkedHashMap<U, Integer> getHashCols() {
//        return hashCols;
//    }
//
//    public void setHashCols(LinkedHashMap<U, Integer> hashCols) {
//        this.hashCols = hashCols;
//    }
//
//    public ArrayList<T> getRowObjects() {
//        return new ArrayList<T>(hashRows.keySet());
//    }
//    
//    public void setRowObjects(ArrayList<T> arrayList) {
//        LinkedHashMap<T, Integer> newHashRows = new LinkedHashMap<T, Integer>((int) Math.ceil(arrayList.size() / 0.75));
//        int i = 0;
//        for(T s: arrayList){
//            if(!newHashRows.containsKey(s)){
//                newHashRows.put(s, i);
//            } else {
//                System.out.println("Error, new row names contains dupilcates.");
//                System.exit(-1);
//            }
//            i++;
//        }
//        
//        this.hashRows = newHashRows;
//    }
//
//    public ArrayList<U> getColObjects() {
//        return new ArrayList<U>(hashCols.keySet());
//    }
//    
//    public void setColObjects(ArrayList<U> arrayList) {
//        LinkedHashMap<U, Integer> newHashCols = new LinkedHashMap<U, Integer>((int) Math.ceil(arrayList.size() / 0.75));
//        int i = 0;
//        for(U s: arrayList){
//            if(!newHashCols.containsKey(s)){
//                newHashCols.put(s, i);
//            } else {
//                System.out.println("Error, new row names contains dupilcates.");
//                System.exit(-1);
//            }
//            i++;
//        }
//         this.hashCols = newHashCols;
//    }
//
//    public int getNrRows() {
//        return nrRows;
//    }
//
//    public void setNrRows(int nrRows) {
//        this.nrRows = nrRows;
//    }
//
//    public int getNrCols() {
//        return nrCols;
//    }
//
//    public void setNrCols(int nrCols) {
//        this.nrCols = nrCols;
//    }
//
//    public DenseDoubleMatrix2D getMatrix() {
//        return Matrix;
//    }
//
//    public void setMatrix(DenseDoubleMatrix2D Matrix) {
//        this.Matrix = Matrix;
//    }
//    
//    public void setMatrix(double[][] Matrix) {
//        this.Matrix = new DenseDoubleMatrix2D(Matrix);
//    }
}
