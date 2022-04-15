package mbqtl.matrix;

import cern.colt.matrix.AbstractMatrix2D;

/**
 * Created by hwestra on 5/16/16.
 */
public class ShortMatrix2D extends AbstractMatrix2D {

	private short[] matrix;

	public ShortMatrix2D(short[][] matrix) {
		rows = matrix.length;
		columns = matrix[0].length;
		this.matrix = new short[rows * columns];
		for (int i = 0; i < rows; i++) {
			int start = i * columns;
			System.arraycopy(matrix[i], 0, this.matrix, start, columns);
		}
	}

	public ShortMatrix2D(int rows, int columns) {
		this.rows = rows;
		this.columns = columns;
		this.matrix = new short[rows * columns];
	}

	public short getQuick(int i, int j) {
		int index = i * columns + j;
		return matrix[index];
	}

	public void setQuick(int i, int j, short b) {
		int index = i * columns + j;
		matrix[index] = b;
	}

	public short[][] toArray() {
		short[][] output = new short[rows][columns];
		for (int i = 0; i < rows; i++) {
			int start = i * columns;
			System.arraycopy(this.matrix, start, output[i], 0, columns);
		}
		return output;
	}


}
