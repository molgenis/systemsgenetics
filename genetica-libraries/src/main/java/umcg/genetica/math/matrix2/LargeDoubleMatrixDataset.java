/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.matrix2;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseLargeDoubleMatrix2D;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.logging.Level;
import java.util.regex.Pattern;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author MarcJan
 */
public class LargeDoubleMatrixDataset<R, C> extends DoubleMatrixDataset<R, C> {

    private DenseLargeDoubleMatrix2D matrix;

	
	
    public LargeDoubleMatrixDataset() {
    }

	public LargeDoubleMatrixDataset(DenseLargeDoubleMatrix2D matrix, LinkedHashMap<R, Integer> hashRows, LinkedHashMap<C, Integer> hashCols) {
		super(hashRows, hashCols);
		this.matrix = matrix;
        
        this.rows = matrix.rows();
        this.columns= matrix.columns();
	}

    public LargeDoubleMatrixDataset(int nrRows, int nrCols) {
        this(nrRows, nrCols, null);
    }

    public LargeDoubleMatrixDataset(int nrRows, int nrCols, Double initialValue) {
        this.rows = nrRows;
        this.columns= nrCols;
        // runtime type of the arrays will be Object[] but they can only contain T and U elements
        this.setHashRows(new LinkedHashMap<R, Integer>((int) Math.ceil(nrRows / 0.75)));
        this.setHashCols(new LinkedHashMap<C, Integer>((int) Math.ceil(nrCols / 0.75)));

        this.matrix = new DenseLargeDoubleMatrix2D(nrRows, nrCols);

        if (initialValue != null) {
            this.matrix.assign(initialValue);
        }
    }

    public LargeDoubleMatrixDataset(String fileName) throws IOException, Exception {
        this(fileName, "\t");
    }

    public LargeDoubleMatrixDataset(String fileName, String ll) throws IOException, Exception {
        loadDoubleData(fileName, ll);
    }

    protected static LargeDoubleMatrixDataset<String, String> loadDoubleData(String fileName, String delimiter) throws IOException, Exception {
        
        LargeDoubleMatrixDataset<String, String> dataset = new LargeDoubleMatrixDataset<String, String>();
                
        Pattern splitPatern = Pattern.compile(delimiter);

        int columnOffset = 1;

        int[] colIndex;
        TextFile in = new TextFile(fileName, TextFile.R);
        String str = in.readLine(); // header
        String[] data = splitPatern.split(str);

        int tmpCols = (data.length - columnOffset);

        dataset.setHashCols(new LinkedHashMap<String, Integer>((int) Math.ceil(tmpCols / 0.75)));

        colIndex = new int[tmpCols];
        for (int s = 0; s < tmpCols; s++) {
            String colName = data[s + columnOffset];
            if(!dataset.getHashCols().containsKey(colName)){
                dataset.getHashCols().put(colName, s);
            } else {
                dataset.LOGGER.warning("Duplicated column name!");
                throw(doubleMatrixDatasetNonUniqueHeaderException);
            }
            colIndex[s] = s + columnOffset;
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

        dataset.setHashRows(new LinkedHashMap<String, Integer>((int) Math.ceil(tmpRows / 0.75)));

        boolean correctData = true;
        while ((str = in.readLine()) != null) {
            data = splitPatern.split(str);

            if (!dataset.getHashRows().containsKey(data[0])) {
                dataset.getHashRows().put( data[0], row);
            } else {
                LOGGER.warning("Duplicated row name!");
                throw(doubleMatrixDatasetNonUniqueHeaderException);
            }

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
        }
        if (!correctData) {
            LOGGER.warning("Your data contains NaN/unparseable values!");
        }
        in.close();
        
        dataset.rows = tmpRows;
        dataset.columns = tmpCols;
        dataset.matrix = new DenseLargeDoubleMatrix2D(tmpRows, tmpCols);
        dataset.matrix.assign(initialMatrix);

        LOGGER.log(Level.INFO, "''{0}'' has been loaded, nrRows: {1} nrCols: {2}", new Object[]{fileName, dataset.rows(), dataset.columns()});
        return(dataset);
    }

//    protected void loadDoubleDataTokenizer(String fileName, String delimiter) throws IOException {
//        
//        Pattern splitPatern = Pattern.compile(delimiter);
//
//        int columnOffset = 1;
//
//        int[] colIndex;
//        TextFile in = new TextFile(fileName, TextFile.R);
//        String str = in.readLine(); // header
//        String[] data = splitPatern.split(str);
//
//        int tmpCols = (data.length - columnOffset);
//
//        this.setHashCols(new LinkedHashMap<C, Integer>((int) Math.ceil(tmpCols / 0.75)));
//
//        colIndex = new int[tmpCols];
//        for (int s = 0; s < tmpCols; s++) {
//            String colName = data[s + columnOffset];
//            if(!this.getHashCols().containsKey((C)colName)){
//                this.getHashCols().put((C)colName, s);
//            } else {
//                this.LOGGER.warning("Duplicated column name!");
//                throw(doubleMatrixDatasetNonUniqueHeaderException);
//            }
//            colIndex[s] = s + columnOffset;
//        }
//
//        int tmpRows = 0;
//
//        while (in.readLine() != null) {
//            tmpRows++;
//        }
//        in.close();
//
//        matrix = new DenseLargeDoubleMatrix2D(tmpRows, tmpCols);
//        this.rows = tmpRows;
//        this.columns = tmpCols;
//        
//        in.open();
//        in.readLine(); // read header
//        int row = 0;
//
//        this.setHashRows(new LinkedHashMap<R, Integer>((int) Math.ceil(this.rows() / 0.75)));
//
//        boolean correctData = true;
//
//        while ((str = in.readLine()) != null) {
//            StringTokenizer st = new StringTokenizer(str, delimiter);
//            int col = 0;
//            while (st.hasMoreTokens()) {
//                if (col != 0) {
//                    double d;
//                    try {
//                        d = Double.parseDouble(st.nextToken());
//                    } catch (NumberFormatException e) {
//                        correctData = false;
//                        d = Double.NaN;
//                    }
//                    matrix.setQuick(row, col, d);
//                } else {
//                    String key = st.nextToken();
//                    if (!this.getHashRows().containsKey((R) key)) {
//                        this.getHashRows().put((R) key, row);
//                    } else {
//                        LOGGER.warning("Duplicated row name!");
//                        throw(doubleMatrixDatasetNonUniqueHeaderException);
//                    }
//                }
//            }
//            row++;
//        }
//        if (!correctData) {
//            LOGGER.warning("Your data contains NaN/unparseable values!");
//        }
//        in.close();
//        LOGGER.log(Level.INFO, "''{0}'' has been loaded, nrRows: {1} nrCols: {2}", new Object[]{fileName, this.rows(), this.columns()});
//    }
	
	@Override
	public LargeDoubleMatrixDataset<C, R> viewDice(){
		return new LargeDoubleMatrixDataset<C, R>((DenseLargeDoubleMatrix2D) matrix.viewDice(), hashCols, hashRows);
	}
    
    @Override
    public DoubleMatrix2D getMatrix() {
        return matrix;
    }

    public DenseLargeDoubleMatrix2D getLargeDenseMatrix() {
        return matrix;
    }

    @Override
    public void setMatrix(double[][] Matrix) {
        this.rows = Matrix.length;
        this.columns=Matrix[0].length;

        matrix = new DenseLargeDoubleMatrix2D(rows(), columns());
        matrix.assign(Matrix);
    }
    
    public void setMatrix(DenseLargeDoubleMatrix2D Matrix) {
        this.matrix = Matrix;
        this.rows = Matrix.rows();
        this.columns=Matrix.columns();
    }

	@Override
	public double[][] elements() {
		return matrix.elements();
	}

	@Override
	public double getQuick(int i, int i1) {
		return getQuick(i, i1);
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
	protected DoubleMatrix1D like1D(int i, int i1, int i2) {
		return like1D(i, i1, i2);
	}

	@Override
	protected DoubleMatrix2D viewSelectionLike(int[] ints, int[] ints1) {
		//implemented as in wrapper double matrix 2d. Only has protected access.
		throw new InternalError(); // should never be called
	}
}
