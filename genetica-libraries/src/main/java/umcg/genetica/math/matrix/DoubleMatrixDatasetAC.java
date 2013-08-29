/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.matrix;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * @author juha
 */
@Deprecated
public abstract class DoubleMatrixDatasetAC<T, U> {

    public enum LoadLabels {

        LOAD_BOTH, LOAD_ROWS, LOAD_COLUMNS, DONT_LOAD
    };

    public Map<T, Integer> hashRows = new HashMap<T, Integer>();
    public Map<U, Integer> hashCols = new HashMap<U, Integer>();
    
    public List<T> rowObjects = null;
    public List<U> colObjects = null;
    public int nrRows;
    public int nrCols;

    public abstract double[] get(int x);

    public abstract double get(int x, int y);

    public abstract void recalculateHashMaps();
}
