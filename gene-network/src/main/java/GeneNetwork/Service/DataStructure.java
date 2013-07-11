/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package GeneNetwork.Service;

import java.util.ArrayList;

/**
 *
 * @author R.S.N. Fehrmann
 */
public class DataStructure {
    
    private ArrayList<SymmetricMatrix<Short> > d_symmetricMatrixList;
    private ArrayList<String> d_pathToSymmetricMatrixList;
    private ArrayList<String> d_identifierOfSymmetricMatrixList;
    
    public DataStructure() {
        d_symmetricMatrixList = new ArrayList<SymmetricMatrix<Short> >();
    }
    
    public void addSymmetricMatrix(SymmetricMatrix<Short> symmetricMatrix, String identifier, String path)
        throws IllegalArgumentException
    {
        if (symmetricMatrix == null)
            throw new IllegalArgumentException("SymmetricMatrix<Short> symmetricMatrix == null");
        
        if (identifier == null)
            throw new IllegalArgumentException("String identifier == null");
        
        if (path == null)
            throw new IllegalArgumentException("String path == null");
        
        d_symmetricMatrixList.add(symmetricMatrix);
        d_pathToSymmetricMatrixList.add(path);
        d_identifierOfSymmetricMatrixList.add(identifier);
    }
    
    public String getIdentifier(int index) 
        throws ArrayIndexOutOfBoundsException
    {
        if (index > d_identifierOfSymmetricMatrixList.size())
            throw new ArrayIndexOutOfBoundsException();
            
        return d_identifierOfSymmetricMatrixList.get(index);
    }
    
    public String getPath(int index) 
        throws ArrayIndexOutOfBoundsException
    {
        if (index > d_pathToSymmetricMatrixList.size())
            throw new ArrayIndexOutOfBoundsException();
        
        return d_pathToSymmetricMatrixList.get(index);        
    }
    
    public SymmetricMatrix<Short> getSymmetricMatrix(int index)
        throws ArrayIndexOutOfBoundsException    
    {
        if (index > d_symmetricMatrixList.size())
            throw new ArrayIndexOutOfBoundsException();
        
        return d_symmetricMatrixList.get(index);
    }
    
    public int size() {
        return d_symmetricMatrixList.size();
    }
}
