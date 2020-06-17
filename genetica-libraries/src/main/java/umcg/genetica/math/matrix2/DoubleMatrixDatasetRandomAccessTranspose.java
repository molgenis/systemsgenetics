package umcg.genetica.math.matrix2;

import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.File;
import java.io.IOException;
import java.util.stream.IntStream;

public class DoubleMatrixDatasetRandomAccessTranspose {


    public static void main(String[] args) {
        try {
            String matrix = "D:\\tmp\\matrixtest\\testmatrix.txt";
            String matrixbinary = "D:\\tmp\\matrixtest\\testmatrix-binary";
            DoubleMatrixDatasetRandomAccessTranspose t = new DoubleMatrixDatasetRandomAccessTranspose();
            t.createMatrix(matrix, 100, 25);
            DoubleMatrixConverter.TextToBinary(matrix, matrixbinary);
//
            String matrixbinarytp = "D:\\tmp\\matrixtest\\testmatrix-binary-transposed";
            String matrixbinarytp2 = "D:\\tmp\\matrixtest\\testmatrix-binary-transposed2";
            String matrixbinarytptxt = "D:\\tmp\\matrixtest\\testmatrix-binary-transposed.txt";
            String matrixbinarytptxt2 = "D:\\tmp\\matrixtest\\testmatrix-binary-transposed2.txt";
//            t.transpose(matrixbinary, matrixbinarytp);
            t.transposeLargeMatrix(matrixbinary, matrixbinarytp2, 3);
//
            DoubleMatrixConverter.BinaryToText(matrixbinarytp, matrixbinarytptxt);
            DoubleMatrixConverter.BinaryToText(matrixbinarytp2, matrixbinarytptxt2);
//
//            String matrixbinarytptp = "D:\\tmp\\matrixtest\\testmatrix-binary-transposed-transposed";
//            t.transpose(matrixbinarytp, matrixbinarytptp);
//            String matrixbinarytptptxt = "D:\\tmp\\matrixtest\\testmatrix-binary-transposed-transposed.txt";
//
            String matrixbinarytp2tp = "D:\\tmp\\matrixtest\\testmatrix-binary-transposed2-transposed";
            t.transposeLargeMatrix(matrixbinarytp2, matrixbinarytp2tp, 20);
            String matrixbinarytp2tptxt = "D:\\tmp\\matrixtest\\testmatrix-binary-transposed2-transposed.txt";
            DoubleMatrixConverter.BinaryToText(matrixbinarytp2tp, matrixbinarytp2tptxt);
//            DoubleMatrixConverter.BinaryToText(matrixbinarytptp, matrixbinarytptptxt);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public void createMatrix(String output, int sizeX, int sizeY) throws IOException {
        double[][] matrix = new double[sizeX][sizeY];
        int ctr = 0;
        String[] rows = new String[sizeX];
        String[] cols = new String[sizeY];

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < sizeY; j++) {
                matrix[i][j] = ctr;
                ctr++;
            }
            rows[i] = "Row-" + i;

        }

        for (int i = 0; i < sizeY; i++) {
            cols[i] = "Col-" + i;
        }

        TextFile tf = new TextFile(output, TextFile.W);
        String header = "-\t" + Strings.concat(cols, Strings.tab);
        tf.writeln(header);
        for (int r = 0; r < sizeX; r++) {
            tf.writeln("Row-" + r + "\t" + Strings.concat(matrix[r], Strings.tab));
        }
        tf.close();

    }

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
        ProgressBar pb = new ProgressBar(reader.nrRows, "Transposing");
        for (int i = 0; i < reader.nrRows; i++) {
            double[] row = reader.getNextRow();
//			System.out.println(Strings.concat(row, Strings.tab));
            for (int j = 0; j < reader.nrCols; j++) {

                writer.write(j, i, row[j]);
            }
            pb.set(i);
        }
        writer.close();
        pb.close();
        reader.close();
    }

    public void transposeLargeMatrix(String in, String out, int rowsToProcessAtOnce) throws IOException {
//        System.out.println("This code is broken");
//        System.exit(-1);

        System.out.println("Transposing big matrix.");
        System.out.println("Input: " + in);
        System.out.println("Output: " + out);


        // Gpio.delete(out);
        DoubleMatrixDatasetRandomAccessReader reader = new DoubleMatrixDatasetRandomAccessReader(in);
        DoubleMatrixDatasetRandomAccessWriter writer = new DoubleMatrixDatasetRandomAccessWriter();
        writer.initializeFullMatrix(reader.getColObjects(), reader.getRowObjects(), out);
        writer.close();
        writer.open(out);


//		System.out.println("Stuff is printed here");
        if (rowsToProcessAtOnce > reader.rows()) {
            rowsToProcessAtOnce = reader.rows();
        }
        double[][] buffer = new double[rowsToProcessAtOnce][];
        double[][] transpose = new double[reader.getColObjects().size()][rowsToProcessAtOnce];
        long nrItems = ((long) reader.getColObjects().size() * reader.getColObjects().size());
        long blocksize = (long) reader.getColObjects().size() * 8 * rowsToProcessAtOnce;
        System.out.println("Block size: " + rowsToProcessAtOnce + " x " + reader.getColObjects().size() + ", " + Gpio.humanizeFileSize(blocksize));
        int ctr = 0;

        int buffernr = 0;
        int nrread = 0;


        ProgressBar pb = new ProgressBar(reader.nrRows, "Reading..");
        for (int i = 0; i < reader.nrRows; i++) {
            if (ctr == rowsToProcessAtOnce) {
                System.out.println();
                // buffer full

                double[][] finalTranspose = transpose;
                IntStream.range(0, buffer.length).parallel().forEach(q -> {
                            for (int z = 0; z < buffer[q].length; z++) {
                                finalTranspose[z][q] = buffer[q][z];
                            }
                        }
                );

                // data has all columns
                int colstart = buffernr * rowsToProcessAtOnce;
                ProgressBar pb2 = new ProgressBar(transpose.length, "Writing..");
                for (int q = 0; q < transpose.length; q++) {
                    writer.writeBlock(q, colstart, finalTranspose[q]);
                    pb2.iterate();
                }
                pb2.close();

                buffernr++;
                ctr = 0;
                pb = new ProgressBar(reader.nrRows, "Reading..");
            }
            if (ctr != rowsToProcessAtOnce) {
                // space left in the buffer.
                double[] row = reader.getNextRow();
                buffer[ctr] = row;
                ctr++;
                nrread++;
                pb.set(nrread);

            }


        }
        pb.close();
        if (ctr > 0) {
            transpose = new double[reader.getColObjects().size()][ctr];

            // buffer full
            for (int q = 0; q < ctr; q++) {
                for (int z = 0; z < buffer[q].length; z++) {

                    transpose[z][q] = buffer[q][z];
                }
            }

            // data has all columns,
            ProgressBar pb2 = new ProgressBar(transpose.length, "Writing..");
            int colstart = buffernr * rowsToProcessAtOnce;
            for (int q = 0; q < transpose.length; q++) {
                writer.writeBlock(q, colstart, transpose[q]);
                pb2.iterate();
            }
            pb2.close();
        }


        writer.close();
        reader.close();

    }


}
