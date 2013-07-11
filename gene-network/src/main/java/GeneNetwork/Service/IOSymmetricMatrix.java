/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package GeneNetwork.Service;

import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;

/**
 *
 * @author R.S.N. Fehrmann
 */
public class IOSymmetricMatrix<T> 
{   
    T d_type;
    int d_nDiagonalElements;
    int d_nElementsToRead;
    DataInputStream d_dataInputStream;
            

    public IOSymmetricMatrix(File file)
        throws IllegalArgumentException, java.io.FileNotFoundException, 
        java.io.EOFException, java.io.IOException
    {
        if (file == null)
            throw new IllegalArgumentException("argument File file == null");

        FileInputStream fileInputStream = new FileInputStream(file);

        d_dataInputStream = new DataInputStream(fileInputStream);
               
        // Read integer from file to identify the number of diagonal elements
        d_nDiagonalElements = d_dataInputStream.readInt();
        
        // Determine the total number of elements to read from file
        d_nElementsToRead = (int) (((long) d_nDiagonalElements * (long)(d_nDiagonalElements + 1)) / 2l);
                
        // Determine total number of elements actually in file
        int nElementsInFile = (int) (file.length() / this.sizeOfType());
        
        // Check if expected number of elements is equal to the actual number of elements in file
        if (d_nElementsToRead != nElementsInFile)
        {
            d_dataInputStream.close();
            throw new java.io.IOException("Expected " + Integer.toString(d_nElementsToRead) + this.nameOfType() + "s with byte size of " + Integer.toString(this.sizeOfType()) + " in file, found " + Integer.toString(nElementsInFile));
        }
    }
    
    public void readBytes(SymmetricMatrix<Byte> symmetricMatrix)
        throws java.io.IOException
    {
        symmetricMatrix = new SymmetricMatrix<Byte>(d_nDiagonalElements);
        
        for (int i = 0; i < d_nElementsToRead; ++i)
            symmetricMatrix.setValueSequentially(i, d_dataInputStream.readByte());
        
        d_dataInputStream.close();
    }

    public void readCharacters(SymmetricMatrix<Character> symmetricMatrix) 
        throws java.io.IOException
    {
        symmetricMatrix = new SymmetricMatrix<Character>(d_nDiagonalElements);
        
        for (int i = 0; i < d_nElementsToRead; ++i)
            symmetricMatrix.setValueSequentially(i, d_dataInputStream.readChar());
        
        d_dataInputStream.close();
    }

    public void readShorts(SymmetricMatrix<Short> symmetricMatrix) 
        throws java.io.IOException
    {
        symmetricMatrix = new SymmetricMatrix<Short>(d_nDiagonalElements);
        
        for (int i = 0; i < d_nElementsToRead; ++i)
            symmetricMatrix.setValueSequentially(i, d_dataInputStream.readShort());
        
        d_dataInputStream.close();
    }
                
    public void readIntegers(SymmetricMatrix<Integer> symmetricMatrix) 
        throws java.io.IOException
    {
        symmetricMatrix = new SymmetricMatrix<Integer>(d_nDiagonalElements);
        
        for (int i = 0; i < d_nElementsToRead; ++i)
            symmetricMatrix.setValueSequentially(i, d_dataInputStream.readInt());
        
        d_dataInputStream.close();
    }

    public void readLongs(SymmetricMatrix<Long> symmetricMatrix) 
        throws java.io.IOException
    {
        symmetricMatrix = new SymmetricMatrix<Long>(d_nDiagonalElements);
        
        for (int i = 0; i < d_nElementsToRead; ++i)
            symmetricMatrix.setValueSequentially(i, d_dataInputStream.readLong());
        
        d_dataInputStream.close();
    }

    public void readFloats(SymmetricMatrix<Float> symmetricMatrix) 
        throws java.io.IOException
    {
        symmetricMatrix = new SymmetricMatrix<Float>(d_nDiagonalElements);
        
        for (int i = 0; i < d_nElementsToRead; ++i)
            symmetricMatrix.setValueSequentially(i, d_dataInputStream.readFloat());
        
        d_dataInputStream.close();
    }

    public void readDoubles(SymmetricMatrix<Double> symmetricMatrix) 
        throws java.io.IOException
    {
        symmetricMatrix = new SymmetricMatrix<Double>(d_nDiagonalElements);
        
        for (int i = 0; i < d_nElementsToRead; ++i)
            symmetricMatrix.setValueSequentially(i, d_dataInputStream.readDouble());
        
        d_dataInputStream.close();
    }

    private int sizeOfType() {
        if (d_type == byte.class)
            return 1;
        else if (d_type == char.class)
            return 2;
        else if (d_type == short.class)
            return 2;        
        else if (d_type == int.class)
            return 4;
        else if (d_type == long.class)
            return 8;
        else if (d_type == float.class)
            return 4;
        else if (d_type == double.class)
            return 8;
        else
            return 0;
    }
    
    private String nameOfType() {
        if (d_type == byte.class)
            return "byte";
        else if (d_type == char.class)
            return "char";
        else if (d_type == short.class)
            return "short";        
        else if (d_type == int.class)
            return "int";
        else if (d_type == long.class)
            return "long";
        else if (d_type == float.class)
            return "float";
        else if (d_type == double.class)
            return "double";
        else
            return null;
    }
}
