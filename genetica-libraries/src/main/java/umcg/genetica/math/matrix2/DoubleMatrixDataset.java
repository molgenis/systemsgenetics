/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.matrix2;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author MarcJan
 */
public abstract class DoubleMatrixDataset<R, C> extends DoubleMatrix2D{

    static final Logger LOGGER = Logger.getLogger(DoubleMatrixDataset.class.getName());
    protected LinkedHashMap<R, Integer> hashRows;
    protected LinkedHashMap<C, Integer> hashCols;
    protected Pattern splitPatern;

    public DoubleMatrixDataset(){
		hashRows = new LinkedHashMap<R, Integer>();
		hashCols = new LinkedHashMap<C, Integer>();
    }

	public DoubleMatrixDataset(LinkedHashMap<R, Integer> hashRows, LinkedHashMap<C, Integer> hashCols) {
		this.hashRows = hashRows;
		this.hashCols = hashCols;
	}

    public abstract DoubleMatrix2D getMatrix();

    public void saveLowMemory(String fileName) throws IOException {
        TextFile out = new TextFile(fileName, TextFile.W);

        ArrayList<String> colObjects = (ArrayList<String>) hashCols.keySet();
        ArrayList<String> rowObjects = (ArrayList<String>) hashRows.keySet();

        out.append('-');
        for (int s = 0; s < columns; s++) {
            out.append('\t');
            out.append(colObjects.get(s));
        }
        out.append('\n');

        for (int r = 0; r < rows; r++) {
            out.append(rowObjects.get(r));
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

        ArrayList<String> colObjects = new ArrayList(hashCols.keySet());
        ArrayList<String> rowObjects = new ArrayList(hashRows.keySet());

        out.append('-');
        for (int s = 0; s < getMatrix().columns(); s++) {
            out.append('\t');
            out.append(colObjects.get(s));
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

    public void setRowObjects(ArrayList<R> arrayList) {
        LinkedHashMap<R, Integer> newHashRows = new LinkedHashMap<R, Integer>((int) Math.ceil(arrayList.size() / 0.75));
        int i = 0;
        for (R s : arrayList) {
            if (!newHashRows.containsKey(s)) {
                newHashRows.put(s, i);
            } else {
                System.out.println("Error, new row names contains dupilcates.");
                System.exit(-1);
            }
            i++;
        }

        this.hashRows = newHashRows;
    }

    public ArrayList<C> getColObjects() {
        return new ArrayList<C>(hashCols.keySet());
    }

    public void setColObjects(ArrayList<C> arrayList) {
        LinkedHashMap<C, Integer> newHashCols = new LinkedHashMap<C, Integer>((int) Math.ceil(arrayList.size() / 0.75));
        int i = 0;
        for (C s : arrayList) {
            if (!newHashCols.containsKey(s)) {
                newHashCols.put(s, i);
            } else {
                System.out.println("Error, new row names contains dupilcates.");
                System.exit(-1);
            }
            i++;
        }
        this.hashCols = newHashCols;
    }

    public void setNrRows(int nrRows) {
        this.rows = nrRows;
    }

    public void setNrCols(int nrCols) {
        this.columns = nrCols;
    }

    public abstract void setMatrix(double[][] Matrix);
}
