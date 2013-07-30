/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.matrix2;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.StringTokenizer;
import java.util.logging.Level;
import java.util.regex.Pattern;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author MarcJan
 */
public class SmallDoubleMatrixDataset<R, C> extends DoubleMatrixDataset<R, C>{

    private DenseDoubleMatrix2D matrix;
    
    public SmallDoubleMatrixDataset(){
        
    }

    public SmallDoubleMatrixDataset(DenseDoubleMatrix2D matrix, LinkedHashMap<R, Integer> hashRows, LinkedHashMap<C, Integer> hashCols) {
		super(hashRows, hashCols);
		this.matrix = matrix;
        
        this.rows = matrix.rows();
        this.columns= matrix.columns();
	}
    
    public SmallDoubleMatrixDataset(int nrRows, int nrCols) {
        this(nrRows, nrCols, null);
    }
    
    public SmallDoubleMatrixDataset(int nrRows, int nrCols, Double initialValue) {
        this.rows = nrRows;
        this.columns= nrCols;
        // runtime type of the arrays will be Object[] but they can only contain T and U elements
        this.setHashRows(new LinkedHashMap<R, Integer>((int) Math.ceil(nrRows / 0.75)));
        this.setHashCols(new LinkedHashMap<C, Integer>((int) Math.ceil(nrCols / 0.75)));
        
        this.matrix = new DenseDoubleMatrix2D(nrRows, nrCols);
        
        if(initialValue!=null){
            this.matrix.assign(initialValue);
        }
    }
    
    public SmallDoubleMatrixDataset(String fileName) throws IOException {
        this(fileName, "\t");
    }

    public SmallDoubleMatrixDataset(String fileName, String ll) throws IOException {
        loadDoubleData(fileName, ll);
    }

    protected void loadDoubleData(String fileName, String delimiter) throws IOException {
        splitPatern = Pattern.compile(delimiter);

        int columnOffset = 1;

        int[] colIndex;
        TextFile in = new TextFile(fileName, TextFile.R);
        String str = in.readLine(); // header
        String[] data = splitPatern.split(str);

        int tmpCols = (data.length - columnOffset);

        this.setHashCols(new LinkedHashMap<C, Integer>((int) Math.ceil(tmpCols / 0.75)));

        colIndex = new int[tmpCols];
        for (int s = 0; s < tmpCols; s++) {
            String colName = data[s + columnOffset];
            if(!this.getHashCols().containsKey((C)colName)){
                this.getHashCols().put((C)colName, s);
            } else {
                this.LOGGER.warning("Duplicated column name!");
                System.exit(0);
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

        this.setHashRows(new LinkedHashMap<R, Integer>((int) Math.ceil(tmpRows / 0.75)));

        boolean correctData = true;
        while ((str = in.readLine()) != null) {
            data = splitPatern.split(str);
            
            if(!this.getHashRows().containsKey((R)data[0])){
                this.getHashRows().put((R) data[0], row);
            } else {
                LOGGER.warning("Duplicated row name!");
                System.exit(0);
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

        matrix = new DenseDoubleMatrix2D(initialMatrix);
        this.rows = tmpRows;
        this.columns = tmpCols;
        
        LOGGER.log(Level.INFO, "''{0}'' has been loaded, nrRows: {1} nrCols: {2}", new Object[]{fileName, this.rows(), this.columns()});
    }

    protected void loadDoubleDataTokenizer(String fileName, String delimiter) throws IOException {
        splitPatern = Pattern.compile(delimiter);

        int columnOffset = 1;

        int[] colIndex;
        TextFile in = new TextFile(fileName, TextFile.R);
        String str = in.readLine(); // header
        String[] data = splitPatern.split(str);

        int tmpCols = (data.length - columnOffset);

        this.setHashCols(new LinkedHashMap<C, Integer>((int) Math.ceil(tmpCols / 0.75)));

        colIndex = new int[tmpCols];
        for (int s = 0; s < tmpCols; s++) {
            String colName = data[s + columnOffset];
            if(!this.getHashCols().containsKey((C)colName)){
                this.getHashCols().put((C)colName, s);
            } else {
                this.LOGGER.warning("Duplicated column name!");
                System.exit(0);
            }
            colIndex[s] = s + columnOffset;
        }

        int tmpRows = 0;

        while (in.readLine() != null) {
            tmpRows++;
        }
        in.close();

        matrix = new DenseDoubleMatrix2D(tmpRows, tmpCols);
        this.rows = tmpRows;
        this.columns = tmpCols;

        in.open();
        in.readLine(); // read header
        int row = 0;

        this.setHashRows(new LinkedHashMap<R, Integer>((int) Math.ceil(this.rows / 0.75)));

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
                    matrix.setQuick(row, col, d);
                } else {
                    String key = st.nextToken();
                    if(!this.getHashRows().containsKey((R) key)){
                        this.getHashRows().put((R)key, row);
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
        LOGGER.log(Level.INFO, "''{0}'' has been loaded, nrRows: {1} nrCols: {2}", new Object[]{fileName, this.rows(), this.columns()});
    }

	@Override
	public SmallDoubleMatrixDataset<C, R> viewDice(){
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
        this.columns=Matrix[0].length;
    }
    
    public void setMatrix(DenseDoubleMatrix2D Matrix) {
        this.matrix = Matrix;
        this.rows = Matrix.rows();
        this.columns=Matrix.columns();
    }
	
	@Override
	public double[] elements() {
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
		//this cannot be implemented
		throw new InternalError();
	}
}
