/*
 * SymmetricShortDistancaMatrix.java
 *
 * Created on 04 June 2004, 16:03
 */

package umcg.genetica.math.matrix;

import umcg.genetica.io.bin.BinaryFile;

import java.io.IOException;

/**
 * Symmetric Short Distance Matrix: A memory efficient matrix capable of performing calculations on all genes x all genes:
 *
 * @author Lude Franke
 */
public class SymmetricFloatDistanceMatrix {

    private int size;
    private float[] matrix;
    private final static double eLog10 = Math.log(10);
    private final static int MAX_ALL_PAIRS_STEPS = 50;
    public final static float DEFAULT_VALUE = Float.NaN;
    private long[] elementIndex;

    /**
     * Creates a new instance of SymmetricShortDistancaMatrix
     * <p>
     * Defines a new Symmetric matrix, with a predefined size
     * <p>
     * Initially all values will be Short.MAX_Value (65535)
     * <p>
     * All data is stored in memory efficient one dimensional array, which costs ((size * (size + 1)) / 2) * 4 bytes (float)
     */
    public SymmetricFloatDistanceMatrix(int size) {

        this.size = size;
        long arraySize = ((long) size * (long) (size + 1)) / 2l;
        matrix = new float[(int) arraySize];
        elementIndex = new long[size];
        for (int x = 0; x < size; x++) elementIndex[x] = (long) x * (long) size - ((long) x * (long) (x + 1)) / 2l;
        setDefaultValue();
    }

    public SymmetricFloatDistanceMatrix(int size, boolean setMaxDistance) {

        this.size = size;
        long arraySize = ((long) size * (long) (size + 1)) / 2l;
        matrix = new float[(int) arraySize];
        elementIndex = new long[size];
        for (int x = 0; x < size; x++) elementIndex[x] = (long) x * (long) size - ((long) x * (long) (x + 1)) / 2l;
        if (setMaxDistance) setDefaultValue();
    }

    public void setAllElements(float value) {
        for (int x = 0; x < size; x++) {
            for (int y = x; y < size; y++) {
                matrix[getElementIndex(x, y)] = value;
            }
        }
    }

    private void setDefaultValue() {
        setAllElements(DEFAULT_VALUE);
    }

    /* Get the array element location for point (x,y):
     */
    private int getElementIndex(int x, int y) {
        if (x > y) {
            return (int) elementIndex[y] + x;
        } else {
            return (int) elementIndex[x] + y;
        }
    }

    /* Set the value for point (x,y):
     */
    public void set(int x, int y, float value) {
        matrix[getElementIndex(x, y)] = value;
    }

    /* Get the value for point (x,y):
     */
    public float get(int x, int y) {
        return matrix[getElementIndex(x, y)];
    }

    /* Get size of matrix:
     */
    public int size() {
        return size;
    }

    /* Get the shortest path for all the pairs using the Floyd-Warshall algorithm:
     *
     * Please be aware that this method replaces the values inside the current matrix with the shortest path values!
     */
    public void save(java.io.File fileName) throws IOException {

        BinaryFile bf = new BinaryFile(fileName.toString(), BinaryFile.W);
        for (int i = 0; i < matrix.length; i++) {
            bf.writeFloat(matrix[i]);
        }
        bf.close();

    }

    public void load(java.io.File fileName) throws IOException {


        long length = fileName.length();
        long matrixLength = 4l * matrix.length; // 4 bytes per element
        if (length == matrixLength) {
            BinaryFile bf = new BinaryFile(fileName.toString(), BinaryFile.R);

            for (int i = 0; i < matrix.length; i++) {
                matrix[i] = bf.readFloat();
            }

            bf.close();
        }
    }
}