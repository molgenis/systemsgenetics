package umcg.genetica.math.matrix2;

import umcg.genetica.text.Strings;

import java.io.File;
import java.io.IOException;

public class DoubleMatrixDatasetRandomAccessTranspose {

	public void transpose(String in, String out) throws IOException {
		System.out.println("Transpose.");
		System.out.println(in);
		System.out.println(out);
		// Gpio.delete(out);
		DoubleMatrixDatasetRandomAccessReader reader = new DoubleMatrixDatasetRandomAccessReader(in);
		DoubleMatrixDatasetRandomAccessWriter writer = new DoubleMatrixDatasetRandomAccessWriter();
		writer.initializeFullMatrix(reader.getColObjects(), reader.getRowObjects(), out);
		writer.close();
		writer.open(out);
		for (int i = 0; i < reader.nrRows; i++) {
			double[] row = reader.getNextRow();
			System.out.println(Strings.concat(row, Strings.tab));
			for (int j = 0; j < reader.nrCols; j++) {

				writer.write(j, i, row[j]);
			}
		}
		writer.close();
	}

	public void transposeLargeMatrix(String in, String out, int rowsToProcessAtOnce) throws IOException {
		System.out.println("Transposing big matrix.");
		System.out.println(in);
		System.out.println(out);

		// Gpio.delete(out);
		DoubleMatrixDatasetRandomAccessReader reader = new DoubleMatrixDatasetRandomAccessReader(in);
		DoubleMatrixDatasetRandomAccessWriter writer = new DoubleMatrixDatasetRandomAccessWriter();
		writer.initializeFullMatrix(reader.getColObjects(), reader.getRowObjects(), out);
		writer.close();
		writer.open(out);


		System.out.println("Stuff is printed here");
		if (rowsToProcessAtOnce > reader.rows()) {
			rowsToProcessAtOnce = reader.rows();
		}
		double[][] buffer = new double[rowsToProcessAtOnce][];
		double[][] transpose = new double[reader.getColObjects().size()][rowsToProcessAtOnce];

		int ctr = 0;

		int buffernr = 0;
		for (int i = 0; i < reader.nrRows; i++) {
			if (ctr == rowsToProcessAtOnce) {
				// buffer full
				for (int q = 0; q < buffer.length; q++) {
					for (int z = 0; z < buffer[q].length; z++) {
						transpose[z][q] = buffer[q][z];
					}
				}

				// data has all columns
				int colstart = buffernr * rowsToProcessAtOnce;
				for (int q = 0; q < transpose.length; q++) {
					writer.writeBlock(q, colstart, transpose[q]);
				}

				buffernr++;
				ctr = 0;
			}
			if (ctr != rowsToProcessAtOnce) {
				// space left in the buffer.
				double[] row = reader.getNextRow();
				buffer[ctr] = row;
				ctr++;

			}


		}

		if (ctr > 0) {
			transpose = new double[reader.getColObjects().size()][ctr];

			// buffer full
			for (int q = 0; q < ctr; q++) {
				for (int z = 0; z < buffer[q].length; z++) {

					transpose[z][q] = buffer[q][z];
				}
			}

			// data has all columns,
			int colstart = buffernr * rowsToProcessAtOnce;
			for (int q = 0; q < transpose.length; q++) {
				writer.writeBlock(q, colstart, transpose[q]);
			}
		}
		writer.close();


	}


}
