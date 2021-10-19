package umcg.genetica.math.matrix2;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.stream.IntStream;

public class DoubleMatrixConverter {

    public static void TextToBinary(String in, String out) throws IOException {
        System.out.println("Input: " + in);
        System.out.println("Output: " + out);
        TextFile tf = new TextFile(in, TextFile.R);
        ArrayList<String> colids = new ArrayList<>();
        String[] header = tf.readLineElems(Strings.whitespace);
        String[] elems = tf.readLineElems(Strings.whitespace);
        int starti = 1;
        int nrcols = header.length - 1;
        if (header.length == elems.length - 1) {
            starti = 1;
            nrcols = header.length;
        }

        for (int i = starti; i < header.length; i++) {
            colids.add(header[i]);
        }

        System.out.println("Matrix has " + colids.size() + " cols.");
        int printevery = 100;
        if (colids.size() > 1000000) {
            printevery = 1;
        }
        DoubleMatrixDatasetAppendableWriter writer = new DoubleMatrixDatasetAppendableWriter(colids, out);
        double[] ln = new double[nrcols];
        int lnctr = 0;

        while (elems != null) {

            if (elems.length != nrcols + 1) {
                System.err.println("Error on line " + lnctr + " with row id: " + elems[0] + ". Expected " + nrcols + " columns, but found " + (elems.length - 1));
                System.exit(-1);
            } else {
                String rowid = elems[0];
                String[] finalElems = elems;
                IntStream.range(1, elems.length).parallel().forEach(i -> {
                    ln[i - 1] = Double.parseDouble(finalElems[i]);
                });

                writer.append(ln, rowid);

            }
            lnctr++;


            if (lnctr % printevery == 0) {
                System.out.print(lnctr + "\tlines processed\r");
            }

            elems = tf.readLineElems(Strings.whitespace);
        }
        System.out.print(lnctr + "\tlines processed. Done\n");
        tf.close();
        writer.close();
    }

    public static void BinaryToText(String in, String out) throws IOException {

        DoubleMatrixDatasetRowIterable it = new DoubleMatrixDatasetRowIterable(in);
        ArrayList<String> colids = new ArrayList<>();
        colids.add("-");
        colids.addAll(it.getCols());
        TextFile tf = new TextFile(out, TextFile.W);
        tf.writeln(Strings.concat(colids, Strings.tab));

        int rctr = 0;
        ArrayList<String> rowids = new ArrayList<>(it.getRows());
        for (double[] row : it) {
            String lnout = rowids.get(rctr) + "\t" + Strings.concat(row, Strings.tab);
            tf.writeln(lnout);
            rctr++;
        }

        it.close();
        tf.close();
    }


}
