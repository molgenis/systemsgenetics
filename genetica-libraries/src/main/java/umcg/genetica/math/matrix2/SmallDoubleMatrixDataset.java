/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.matrix2;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import java.util.LinkedHashMap;

/**
 *
 * @author MarcJan
 */
public class SmallDoubleMatrixDataset<R, C> extends DoubleMatrixDataset<R, C> {

    private DenseDoubleMatrix2D matrix;

    public SmallDoubleMatrixDataset() {
    }

    public SmallDoubleMatrixDataset(DenseDoubleMatrix2D matrix, LinkedHashMap<R, Integer> hashRows, LinkedHashMap<C, Integer> hashCols) {
        super(hashRows, hashCols);
        this.matrix = matrix;

        this.rows = matrix.rows();
        this.columns = matrix.columns();
    }

    public SmallDoubleMatrixDataset(int nrRows, int nrCols) {
        this(nrRows, nrCols, null);
    }

    public SmallDoubleMatrixDataset(int nrRows, int nrCols, Double initialValue) {
        this.rows = nrRows;
        this.columns = nrCols;
        // runtime type of the arrays will be Object[] but they can only contain T and U elements
        this.setHashRows(new LinkedHashMap<R, Integer>((int) Math.ceil(nrRows / 0.75)));
        this.setHashCols(new LinkedHashMap<C, Integer>((int) Math.ceil(nrCols / 0.75)));

        this.matrix = new DenseDoubleMatrix2D(nrRows, nrCols);

        if (initialValue != null) {
            this.matrix.assign(initialValue);
        }
    }

    @Override
    public SmallDoubleMatrixDataset<C, R> viewDice() {
        return new SmallDoubleMatrixDataset<C, R>((DenseDoubleMatrix2D) matrix.viewDice(), hashCols, hashRows);
    }

    @Override
    public DoubleMatrix2D getMatrix() {
        return matrix;
    }

    public DenseDoubleMatrix2D getSmallDenseMatrix() {
        return matrix;
    }

    @Override
    public void setMatrix(double[][] Matrix) {
        matrix = new DenseDoubleMatrix2D(Matrix);
        this.rows = Matrix.length;
        this.columns = Matrix[0].length;
    }

    public void setMatrix(DenseDoubleMatrix2D Matrix) {
        this.matrix = Matrix;
        this.rows = Matrix.rows();
        this.columns = Matrix.columns();
    }

    @Override
    public double[] elements() {
        return matrix.elements();
    }

    @Override
    public double getQuick(int i, int i1) {
        return matrix.getQuick(i, i1);
    }

    @Override
    public DoubleMatrix2D like(int i, int i1) {
        return matrix.like(i, i1);
    }

    @Override
    public DoubleMatrix1D like1D(int i) {
        return matrix.like1D(i);
    }

    @Override
    public void setQuick(int i, int i1, double d) {
        matrix.setQuick(i, i1, d);
    }

    @Override
    public DoubleMatrix1D vectorize() {
        return matrix.vectorize();
    }

    @Override
    protected DoubleMatrix2D viewSelectionLike(int[] ints, int[] ints1) {
        //this cannot be implemented
        throw new InternalError();
    }
}
