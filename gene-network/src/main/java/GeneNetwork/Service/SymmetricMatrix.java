/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package GeneNetwork.Service;

import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * @author R.S.N. Fehrmann
 */
public class SymmetricMatrix<T> {
    
    private int d_numberOfDiagonalElements;
    private ArrayList<T> d_matrix;
    private int[] d_translateTable;
    private HashMap<String, Integer> d_identifierToIndex;
    
    public SymmetricMatrix() {    
    }
    
    public SymmetricMatrix(int numberOfDiagonalElements) {
        initMatrix(numberOfDiagonalElements);
    }

    public SymmetricMatrix(ArrayList<String> identifiers)
    {
        for (int i = 0; i < identifiers.size(); ++i) 
            d_identifierToIndex.put(identifiers.get(i), i);
    
        initMatrix(identifiers.size());
    }
    
    public void initMatrix(int n) {

        d_numberOfDiagonalElements = n;
  
        d_matrix = new ArrayList<T>((int) ((long) n * (long)(n + 1) / 2l));
  
        d_translateTable = new int[n];

        for (int x = 0; x < n; ++x) 
            d_translateTable[x] = (int) ((long) x * (long) n - ((long) x * (long)(x + 1)) / 2l);
    }
    
    T getValue(int x, int y) {
        if (x > y) {
            return d_matrix.get(d_translateTable[y] + x);
        } else {
            return d_matrix.get(d_translateTable[x] + y);
        }    
    }

    T getValue(String strA, String strB) {
        return getValue(d_identifierToIndex.get(strA), d_identifierToIndex.get(strB));
    }
    
    void setValue(int x, int y, T value) {
        if (x > y) {
            d_matrix.set(d_translateTable[y] + x, value);
        } else {
            d_matrix.set(d_translateTable[x] + y, value);
        }    
    }

    void setValue(String strA, String strB, T value) {
        setValue(d_identifierToIndex.get(strA), d_identifierToIndex.get(strB), value);
    }

    void setValueSequentially(int index, T value) {
        d_matrix.set(index, value);
    }
    
    boolean containsIdentifier(String identifier) {
        return d_identifierToIndex.containsKey(identifier);
    }
    
    public int getNumberOfDiagonalElements() {
        return d_numberOfDiagonalElements;
    }
    
    public void clear() {
        d_matrix.clear();
        d_translateTable = null;
    }
}
